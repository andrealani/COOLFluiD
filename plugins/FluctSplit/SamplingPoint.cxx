#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/SamplingPoint.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SamplingPoint,
                      DataProcessingData,
		      FluctSplitSpaceTimeModule>
aSamplingPointProvider("SamplingPoint");

//////////////////////////////////////////////////////////////////////////////

void SamplingPoint::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("OutputFile","Name of Output File.");
  options.addConfigOption< std::vector<CFreal> >("location","Location of the probe.");
  options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file.");
}

//////////////////////////////////////////////////////////////////////////////

SamplingPoint::SamplingPoint( const std::string& name) :
  DataProcessingCom(name),
  m_fhandle(),
  socket_nodes("nodes"),
  socket_states("states"),
  m_location(0)
{
  addConfigOptionsTo(this);

  m_OutputFile = "Probe.dat";
  setParameter("OutputFile",&m_OutputFile);

  setParameter("location" ,&m_location);

  m_SaveRate = 1;
  setParameter("SaveRate",&m_SaveRate);
}

//////////////////////////////////////////////////////////////////////////////

SamplingPoint::~SamplingPoint()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SamplingPoint::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SamplingPoint::setup()
{
  CFAUTOTRACE;

  DataProcessingCom::setup(); // first call setup of parent class

  m_cellbuilder.setup(); // Set up this thing

  RealVector coords ( m_location.size() );

  coords[XX] = m_location[XX];
  coords[YY] = m_location[YY];

  m_success = false;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner TRS's
  std::vector<std::string> tags; tags.push_back ("inner"); tags.push_back ("cell");
  std::vector< Common::SafePtr<TopologicalRegionSet> > innertrs =
     MeshDataStack::getActive()->getFilteredTrsList(tags);

  // loop over the inner TRS's
  StdTrsGeoBuilder::GeoData& geoData = m_cellbuilder.getDataGE();
  for (CFuint iTRS = 0; iTRS < innertrs.size(); ++iTRS)
  {
    //SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");

    //prepares to loop over cells by getting the GeometricEntityPool
    geoData.trs = innertrs[iTRS];//trs;

    // loop over element types
    for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
    {
      // get the number of elements
      const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

      // get start index of this element type in global element list
      CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

      // loop over elements
      for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
      {
        // build the GeometricEntity
        geoData.idx = cellIdx;

        GeometricEntity *const cell = m_cellbuilder.buildGE();

        // serach algorithm to find the cell which contains the point
        bool found = cell->isInElement ( coords );

        if ( found )
        {
           m_success = true;
           m_SamplingCellIndex = cellIdx;
           m_PointerToSampledTRS = innertrs[iTRS];//trs;  // store the cellIdx, trs pointer ( cellIdx, trs )
        }

        m_cellbuilder.releaseGE(); // always release the cell before going to the next one

      } //end Loop over elements

    }//end Loop over element types

  } // end loop TRS's


  //////////////////

  if (m_success)
  {
    m_fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

    cf_assert (!m_fhandle->isopen());

    boost::filesystem::path m_path_name = PathAppender::getInstance().appendParallel(m_OutputFile);
    ofstream& myfile = m_fhandle->open(m_path_name);

    myfile << "TITLE = \"Probe values\"" << "\n";
    myfile << "FILETYPE = FULL" << "\n";
    myfile << "VARIABLES = \"t\" \"x0\" \"x1\" \"state0\" \"state1\" \"state2\" \"state3\"" << "\n";
    myfile << "ZONE T = " << " ";
    myfile << m_OutputFile << "\n";
    myfile << "DATAPACKING = POINT" << "\n";

    m_fhandle->close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SamplingPoint::execute()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();

  const CFuint iter = ssys_status->getNbIter();

  if ( !(iter % m_SaveRate) && (m_success)) { TakeSample(); }
}

//////////////////////////////////////////////////////////////////////////////

void SamplingPoint::TakeSample()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFreal time = subSysStatus->getCurrentTime();

  //RealVector coords ( m_location.size() );
  //coords[XX] = m_location[XX];
  //coords[YY] = m_location[YY];

//------------------------------

  StdTrsGeoBuilder::GeoData& geoData = m_cellbuilder.getDataGE();

  geoData.trs = m_PointerToSampledTRS;//m_found_trs;
  geoData.idx = m_SamplingCellIndex;//m_found_cell_idx;

  GeometricEntity *const cell = m_cellbuilder.buildGE();

  std::vector<State*>* cellStates = cell->getStates();

  const CFuint nbStates = cellStates->size();

  RealVector sf (nbStates);

  RealVector coords( m_location.size() );
  coords[XX] = m_location[XX];
  coords[YY] = m_location[YY];

  const RealVector xsi = cell->computeMappedCoordFromCoord( coords );
  /// up until here the code can be executed in setup()

  sf = cell->computeShapeFunctionAtMappedCoord( xsi );

  // Variable that stores the interpolated state
  State interpolated_state;

  for(CFuint istate = 0; istate < nbStates; ++istate)
  {
    const State& lstate = *(*cellStates)[istate];
    interpolated_state += sf[istate] * lstate;
  }

  m_cellbuilder.releaseGE(); // always release the cell

  const CFreal x_p = coords[XX];
  const CFreal y_p = coords[YY];
  const CFreal rho = interpolated_state[0];
  const CFreal rho0u = interpolated_state[1];
  const CFreal rho0v = interpolated_state[2];
  const CFreal p = interpolated_state[3];

//------------------------------
  m_fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  boost::filesystem::path m_path_name = PathAppender::getInstance().appendParallel(m_OutputFile);

  cf_assert ( !m_fhandle->isopen() );

  ofstream& myfile = m_fhandle->open(m_path_name, ios::app);

  myfile
  << time << " "
  << x_p << " "
  << y_p << " "
  << rho  << " "
  << rho0u << " "
  << rho0v << " "
  << p    << "\n";

  m_fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void SamplingPoint::unsetup()
{
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void SamplingPoint::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );

  if ( ! ( m_location.size() == DIM_2D ||  m_location.size() == DIM_3D ) )
    throw Common::BadValueException ( FromHere(), "SmaplingPoint coordinates must have 2 or 3 dimensions\n" );
}

//////////////////////////////////////////////////////////////////////////////

    } //namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

