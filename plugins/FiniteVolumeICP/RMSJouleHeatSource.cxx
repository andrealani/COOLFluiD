#include "Common/PE.hh"
#include "MathTools/MathConsts.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/BadValueException.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeICP/RMSJouleHeatSource.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RMSJouleHeatSource, DataProcessingData, FiniteVolumeICPModule>
rmsJouleHeatSourceProvider("RMSJouleHeatSource");

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Freq","Frequency frequency of the torch [MHz].");
   options.addConfigOption< CFreal,Config::DynamicOption<> >("DesiredPower", "Desired torch power [kW].");
   options.addConfigOption< CFreal >("Permeability","Permeability of the free space.");
   options.addConfigOption< CFuint >("NbCoils","Number of coils.");
   options.addConfigOption< std::vector<CFreal> >("RadiusCoils","Radius of each circular coil.");
   options.addConfigOption< std::vector<CFreal> >("ZPositionCoils","Position on z-axis of each coil.");
   options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
   options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix _time#.");
   options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix _iter#.");
   options.addConfigOption< std::string >("OutputFileElCurrent","Name of Output File to write the electric current.");
   options.addConfigOption< std::string >("OutputFileEMField","Name of Output File to write the electromagnetic field.");
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSource::RMSJouleHeatSource(const std::string& name) :
  DataProcessingCom(name),
  socket_rmsJouleHeatSource("rmsJouleHeatSource"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  m_library(CFNULL),
  m_geoBuilder(),
  m_physicalData(),
  m_rescaleElCurrent(0.),
  m_elCurrent(1.),
  m_elFieldReInNodes(CFNULL),
  m_elFieldImInNodes(CFNULL),
  m_rmsJouleHeatSourceInNodes(CFNULL)
{
  addConfigOptionsTo(this);
  
  m_freqMHz = 0.;
  setParameter("Freq",&m_freqMHz);
  
  m_desiredPowerkW = 0.;
  setParameter("DesiredPower",&m_desiredPowerkW);
  
  m_permeability = 0.;
  setParameter("Permeability",&m_permeability);
  
  m_nbCoils = 0;
  setParameter("NbCoils",&m_nbCoils);
  
  m_radiusCoils = std::vector<CFreal>();
  setParameter("RadiusCoils",&m_radiusCoils);
  
  m_zPositionCoils = std::vector<CFreal>();
  setParameter("ZPositionCoils",&m_zPositionCoils);
  
  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);
  
  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);
  
  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
  
  m_nameOutputFileElCurrent = "ElectricCurrent.plt";
  setParameter("OutputFileElCurrent",&m_nameOutputFileElCurrent);
  
  m_nameOutputFileEMField = "EMField.plt";
  setParameter("OutputFileEMField",&m_nameOutputFileEMField);
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSource::~RMSJouleHeatSource()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
RMSJouleHeatSource::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_rmsJouleHeatSource);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
RMSJouleHeatSource::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::setup()
{
  CFAUTOTRACE;

  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();

  DataHandle<CFreal> rmsJouleHeatSource = socket_rmsJouleHeatSource.getDataHandle();
  rmsJouleHeatSource.resize(nbCells);
  rmsJouleHeatSource = 0.0;

  m_geoBuilder.setup();

  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(m_library.isNotNull());
  
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_physicalData);
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  if (m_radiusCoils.size() != m_nbCoils) {
    throw BadValueException (FromHere(),"Number of coils doesn't match the size of the vector of coils radius");
  }

  if (m_zPositionCoils.size() != m_nbCoils) {
    throw BadValueException (FromHere(),"Number of coils doesn't match the size of the vector of z-axis positions of coils");
  }

  if (m_desiredPowerkW < 1.0e-8) {
      throw BadValueException (FromHere(),"The desired ICP power is zero or negative");
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal RMSJouleHeatSource::ellipticIntegralFirstKind(CFreal const& k)
{
  CFreal elIntFirstKind = 0.;
  const CFreal P = 1.- k*k;

  if(P < 1.e-15) elIntFirstKind = exp(1000.);
  else {
    elIntFirstKind = 1.38629436 + P*(0.096663443 + P*(0.035900924
         + P*(0.037425637 + 0.014511962*P)))
         - log(P)*(0.5 + P*(0.12498594 + P*(0.068802486
         + P*(0.033283553 + 0.0044178701*P))));
  }

  return elIntFirstKind;
}

//////////////////////////////////////////////////////////////////////////////

CFreal RMSJouleHeatSource::ellipticIntegralSecondKind(CFreal const& k)
{
  CFreal elIntSecondKind = 0.;
  const CFreal P = 1.- k*k;

  if(P < 1.e-15) elIntSecondKind = 1.;
  else {
    elIntSecondKind = 1. + P*(0.44325141 + P*(0.062606012 + P*(0.047573836
         + 0.017365065*P)))
         - log(P)*P*(0.24998368 + P*(0.0920018
         + P*(0.040696975 + 0.0052644964*P)));
  }

  return elIntSecondKind;
}

//////////////////////////////////////////////////////////////////////////////

CFreal RMSJouleHeatSource::ellipticIntegralCombined(CFreal const& k)
{
  CFreal elIntComb = 0.;

  if(k < 1.e-15) elIntComb = 0.;
  else elIntComb = ((2./k) - k)*ellipticIntegralFirstKind(k)
                 - (2./k)*ellipticIntegralSecondKind(k);

  return elIntComb;
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::execute()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  CFLogDebugMin( "RMSJouleHeatSource::execute()" << "\n");

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> rmsJouleHeatSource = socket_rmsJouleHeatSource.getDataHandle();
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  StdTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  cf_assert(updateVarSet.isNotNull());
  SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<EulerTerm>();
  cf_assert(eulerTerm.isNotNull());
  const RealVector& refData = eulerTerm->getReferencePhysicalData();
  
  CFreal k = 0.;
  CFreal r = 0.;
  CFreal z = 0.;
  CFreal rCoil = 0.;
  CFreal zCoil = 0.;
  CFreal vectorPotentialRe = 0.;
  const CFreal vectorPotentialIm = 0.;
  CFreal actualPower = 0.;
  CFreal elFieldRe = 0.;
  CFreal elFieldIm = 0.;
  CFreal elConductivity = 0.;
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbTs = m_library->getNbTempVib() + m_library->getNbTe();
  const CFuint TID = nbEqs-nbTs-2-1;
  const CFuint nbNodes = nodes.size();
  m_elFieldReInNodes.resize(nbNodes);
  m_elFieldReInNodes = 0.;
  m_elFieldImInNodes.resize(nbNodes);
  m_elFieldImInNodes = 0.;
  m_rmsJouleHeatSourceInNodes.resize(nbNodes);
  m_rmsJouleHeatSourceInNodes = 0.;
  RealVector volumeSum(0.,nbNodes);
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();

    // set coordinates of cell centroid
    Node& coordinate = currCell->getState(0)->getCoordinates();
    z = coordinate[XX];
    r = coordinate[YY];

    // computation of the real component of vector potential
    vectorPotentialRe = 0.;
    for (CFuint iCoil = 0; iCoil < m_nbCoils; ++iCoil) {
      rCoil = m_radiusCoils[iCoil];
      zCoil = m_zPositionCoils[iCoil];
      k = sqrt(4.*rCoil*r/((rCoil+r)*(rCoil+r)+(z-zCoil)*(z-zCoil)));
      vectorPotentialRe += sqrt(rCoil/r)*ellipticIntegralCombined(k);
    }
    vectorPotentialRe *= .5*m_permeability*m_elCurrent/MathTools::MathConsts::CFrealPi();
    
    elFieldRe =  2.*MathTools::MathConsts::CFrealPi()*m_freqMHz*1000000.*vectorPotentialIm;
    elFieldIm = -2.*MathTools::MathConsts::CFrealPi()*m_freqMHz*1000000.*vectorPotentialRe;
    
    // set temperature and pressure
    State *currState = currCell->getState(0);
    updateVarSet->computePhysicalData(*currState, m_physicalData);
    CFreal pdim = eulerTerm->getPressureFromState(m_physicalData[EulerTerm::P])*refData[EulerTerm::P];
    CFreal Tdim = m_physicalData[EulerTerm::T]*refData[EulerTerm::T];
    CFreal* tVec = (nbTs > 0) ? &(*currState)[TID] : CFNULL;
    elConductivity = m_library->sigma(Tdim, pdim, tVec);
    
    if(isnan(elConductivity)) { ///test
      CF_DEBUG_OBJ(Tdim);
      CF_DEBUG_OBJ(pdim);
      CF_DEBUG_OBJ(elConductivity);
      abort();
    }
    ///    elConductivity = std::max(0.,((Tdim-6000.)/4000.)*3000.); //approximation

     // time-averaged Joule heat source not scaled to rescaled el. current yet
    rmsJouleHeatSource[iCell] = .5*elConductivity*
      (elFieldRe*elFieldRe + elFieldIm*elFieldIm);

    // take into account only contribution from updatable states
    if (currCell->getState(0)->isParUpdatable()) {
      actualPower += rmsJouleHeatSource[iCell]*volumes[iCell]*2.*MathTools::MathConsts::CFrealPi()*r;
    }

    if(!(iter % m_saveRate)) {
      // used to print electromagnetic field and Joule heat source in nodes
      std::vector<Node*> localNodes = *currCell->getNodes();
      for (CFuint iLocalNode = 0; iLocalNode < localNodes.size(); ++iLocalNode) {
        CFuint nodeID = localNodes[iLocalNode]->getLocalID();
        volumeSum[nodeID] += volumes[iCell];
        m_elFieldReInNodes[nodeID] += elFieldRe*volumes[iCell];
        m_elFieldImInNodes[nodeID] += elFieldIm*volumes[iCell];
        m_rmsJouleHeatSourceInNodes[nodeID] += rmsJouleHeatSource[iCell]*volumes[iCell];
      }
    }

    // Joule heat source multiplied by radius - axisymmetric equations
    rmsJouleHeatSource[iCell] *= r;

    m_geoBuilder.releaseGE();
  }
  
  const std::string nsp = this->getMethodData().getNamespace(); 
  /// @TODO AL: temporary solution ....
  CFdouble totalPower = 0.0;
  if (PE::GetPE().GetProcessorCount(nsp) > 1) {
#ifdef CF_HAVE_MPI
    MPI_Allreduce(&actualPower, &totalPower, 1, MPI_DOUBLE, MPI_SUM,
		  PE::GetPE().GetCommunicator(nsp));
#endif
  }
  else {
    totalPower = actualPower;
  }
  
  // parameter rescaling the electric current
  m_rescaleElCurrent = sqrt(m_desiredPowerkW*1000./totalPower);

  // rescaling Joule heat source to achieve desired power
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    rmsJouleHeatSource[iCell] *= m_rescaleElCurrent*m_rescaleElCurrent;
  }



  // writes electric current
  printElCurrentToFile();

  if(!(iter % m_saveRate)) {
    // to average the electromagnetic field and Joule heat source in nodes by
    // number of neighbouring cells
    for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
      m_elFieldReInNodes[iNode] *= m_rescaleElCurrent/volumeSum[iNode];
      m_elFieldImInNodes[iNode] *= m_rescaleElCurrent/volumeSum[iNode];
      m_rmsJouleHeatSourceInNodes[iNode] *= (m_rescaleElCurrent*m_rescaleElCurrent)/volumeSum[iNode];
    }

    // writes electromagnetic field Joule heat source
    writeOutputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path RMSJouleHeatSource::constructFilename()
{
  const bool isParallel = PE::GetPE().IsParallel ();

  if (isParallel) {
    std::ostringstream fname;
    fname << boost::filesystem::basename(boost::filesystem::path(m_nameOutputFileEMField))
          << "-" << PE::GetPE().GetRank("Default")
          << boost::filesystem::extension(boost::filesystem::path(m_nameOutputFileEMField));
  }
  
  CFLog(VERBOSE, "RMSJouleHeatSource::constructFilename() => Writing Electromagnetic Field to : " << m_nameOutputFileEMField << "\n");
  return boost::filesystem::path(m_nameOutputFileEMField);
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  Electromagnetic Field\n";
  outputFile << "VARIABLES = X Y elFieldRe elFieldIm elFieldValue RMSJouleHeatSource\n";
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::writeOutputFile()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  StdTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint nbCells = cells->getLocalNbGeoEnts();

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(constructFilename());

  prepareOutputFile(outputFile);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();

  outputFile << "ZONE N=" << nbNodes
             << ", E=" << nbCells
             << ", F=FEPOINT"
             << ", ET=" << "QUADRILATERAL"
             << "\n";

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {

    const Node& currNode = *nodes[iNode];

    // Output to File
    outputFile.precision(12);
    outputFile << currNode[XX]
               << " "
               << currNode[YY]
               << " "
               << m_elFieldReInNodes[iNode]
               << " "
               << m_elFieldImInNodes[iNode]
               << " "
               << sqrt(m_elFieldReInNodes[iNode]*m_elFieldReInNodes[iNode]
                       + m_elFieldImInNodes[iNode]*m_elFieldImInNodes[iNode])
               << " "
               << m_rmsJouleHeatSourceInNodes[iNode]
               << "\n";
  }

  // writes connectivity
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();

    std::vector<Node*> nodes = *currCell->getNodes();

    // connectivity
    if (nodes.size() == 4) {
      outputFile << nodes[0]->getLocalID()+1 << " " << nodes[1]->getLocalID()+1 <<
      " " << nodes[2]->getLocalID()+1 << " " << nodes[3]->getLocalID()+1 <<"\n";
    }
    else {
      outputFile << nodes[0]->getLocalID()+1 << " " << nodes[1]->getLocalID()+1 <<
      " " << nodes[2]->getLocalID()+1 << " " << nodes[0]->getLocalID()+1 <<"\n";
    }

    m_geoBuilder.releaseGE();
  }

  CFLog(VERBOSE, "RMSJouleHeatSource::writeOutputFile() => Writing of Electromagnetic Field finished." << "\n");
  outputFile.close();
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::printElCurrentToFile()
{
  CFAUTOTRACE;

  using namespace std;

  boost::filesystem::path fpath (m_nameOutputFileElCurrent);
  fpath = PathAppender::getInstance().appendAllInfo( fpath, m_appendIter, m_appendTime);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(fpath,ios::app);

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFreal actualElCurrent = m_rescaleElCurrent*m_elCurrent;

  if (iter == 1) {
    outputFile << "TITLE  =  Electric Current\n";
    outputFile << "VARIABLES = Iter I\n";
  }
  outputFile << iter << " " << actualElCurrent << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSource::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
