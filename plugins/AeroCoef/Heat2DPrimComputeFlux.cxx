#include "Common/PE.hh"
#include "Common/CFMultiMap.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealMatrix.hh"
#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/PhysicalModel.hh"

#include "AeroCoef/DataProcessingHeat.hh"
#include "AeroCoef/Heat2DPrimComputeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Heat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Heat2DPrimComputeFlux,
                      DataProcessingData,
                      DataProcessingHeatModule>
aHeat2DPrimComputeFluxProvider("Heat2DPrimComputeFlux");

//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("OutputFile","Name of Output File to write the wall values.");
   options.addConfigOption< CFuint >("SaveRate","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");

}

//////////////////////////////////////////////////////////////////////////////

Heat2DPrimComputeFlux::Heat2DPrimComputeFlux(const std::string& name) :
  DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_faceNeighCell("faceNeighCell")
{
  m_file = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_nameOutputFile = "Wall";
  setParameter("OutputFile",&m_nameOutputFile);

  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

}

//////////////////////////////////////////////////////////////////////////////

Heat2DPrimComputeFlux::~Heat2DPrimComputeFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Heat2DPrimComputeFlux::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_faceNeighCell);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::setup()
{
  _model = Framework::PhysicalModelStack::getActive()->getImplementor().d_castTo<HeatPhysicalModel>();

  _gradientsT.resize(PhysicalModelStack::getActive()->getDim());
}
//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::executeOnTrs()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute wall distribution
  if(!(iter % m_saveRate)) { computeWall(); }

}

//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFile(); // file is opened here

  Common::SafePtr<std::vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  cf_assert(m_file->isopen());
  ofstream& fout = m_file->get();

  // builder for standard TRS GeometricEntity's taht will be used for the faces
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderFace;
  geoBuilderFace.setup();
  StdTrsGeoBuilder::GeoData& geoDataFace = geoBuilderFace.getDataGE();
  geoDataFace.trs = getCurrentTRS();
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();

  // builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();

  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  RealVector normal(0.0, DIM_2D);

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace){
    geoDataFace.idx = iFace;
    GeometricEntity & currFace = *geoBuilderFace.buildGE();
    CFuint FaceID = currFace.getID();
    // get the face normal
    std::vector<State*>& states = *(currFace.getStates());
    cf_assert(states.size() >= 2);

    // handle to the neighbor cell
    DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
      faceNeighCell = socket_faceNeighCell.getDataHandle();

    const CFuint cellID = faceNeighCell[FaceID].first;
    const CFuint iFaceLocal = faceNeighCell[FaceID].second;

    geoDataCell.trs = faceNeighCell[FaceID].third;
    geoDataCell.idx = cellID;
    GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

    std::vector<RealVector> m_coords(1);
    m_coords[0].resize(2);

    vector<Framework::State*> * m_states= neighborCell.getStates();
    const CFuint nbStatesInCell = m_states->size();

    m_coords[0] = states[0]->getCoordinates() + states[1]->getCoordinates();
    m_coords[0] *= 0.5;

    normal = neighborCell.computeAvgFaceNormals()[iFaceLocal];
    normal.normalize();

    // Compute the shape function gradient at the projected point
    const std::vector<RealMatrix> cellGradients =
      neighborCell.computeSolutionShapeFunctionGradients(m_coords);

    /// Compute the gradient of T
    _gradientsT = 0.;
    for (CFuint iNode = 0; iNode < nbStatesInCell ; ++iNode){
      for (CFuint iDim =0; iDim < _gradientsT.size() ; ++iDim){
        (_gradientsT)[iDim] += (cellGradients[0])(iNode,iDim) * (*(neighborCell.getState(iNode)))[0] ;
      }
    }

    /// Project the gradient on the face normal
    CFreal flux = 0.;
    for (CFuint iDim = 0; iDim < _gradientsT.size() ; ++iDim){
      flux += _gradientsT[iDim]*normal[iDim];
    }

    RealVector shapeFunctions = currFace.computeShapeFunctionAtCoord(m_coords[0]);
    const CFuint nbNodes = shapeFunctions.size();
    RealVector temperature(1);
    temperature = 0.;
    for (CFuint k = 0; k < nbNodes; ++k)
    {
      temperature[0] += shapeFunctions[k] * ((*((states)[k]))[0]);
    }

    _model->setCurrentZone(faceNeighCell[FaceID].third->getName());
    flux *= _model->getConductivity(m_coords[0], temperature);

    CFreal radiativeFlux = 5.6703*0.00000001*temperature[0]*temperature[0]*temperature[0]*temperature[0];

    // output to File
    fout
    << (m_coords[0])[XX]   << " "
    << (m_coords[0])[YY]   << " "
    << temperature[0] << " "
    << flux   << " "
    << radiativeFlux  << "\n";

    geoBuilderFace.releaseGE();
    geoBuilderCell.releaseGE();
  }

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void Heat2DPrimComputeFlux::prepareOutputFile()
{
  using boost::filesystem::path;

  cf_assert (!m_file->isopen());

  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFile);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );

  ofstream& fout = m_file->open(file);

  fout << "TITLE  =  Values at the Wall" << "\n";
  fout << "VARIABLES = x y T Q Qrad_eps1" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

