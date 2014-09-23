#include "LESDataProcessing.hh"
#include "GradientComputerFVMCC.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "MathTools/MathFunctions.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GradientComputerFVMCC, 
                                  LESProcessingData, 
                                  GradientComputer, 
                                  LESDataProcessingModule> 
gradientComputerFVMCCProvider("GradientComputerFVMCC");

// //////////////////////////////////////////////////////////////////////////////

void GradientComputerFVMCC::defineConfigOptions(Config::OptionList& options)
{
}
		
//////////////////////////////////////////////////////////////////////////////
		
GradientComputerFVMCC::GradientComputerFVMCC(const std::string& name) :
  GradientComputer(name),
  m_primState(),
  m_avState(),
  m_normal()
{
}

//////////////////////////////////////////////////////////////////////////////

void GradientComputerFVMCC::configure ( Config::ConfigArgs& args )
{
  GradientComputer::configure(args);

  m_globalSockets.createSocketSink<Framework::State*>("states");
  m_globalSockets.createSocketSink<Framework::Node*>("nodes");
  m_sockets.createSocketSink<Framework::State*>("gstates");
  m_sockets.createSocketSink<RealVector>("nstates");
  
  m_sockets.createSocketSink<CFreal>("volumes");
  m_sockets.createSocketSink<CFint>("isOutward");
  m_sockets.createSocketSink<CFreal>("normals");
}

//////////////////////////////////////////////////////////////////////////////


GradientComputerFVMCC::~GradientComputerFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////


void GradientComputerFVMCC::setup()
{
  CFAUTOTRACE;
  CFLog(INFO, " +++ GradientComputerFVMCC::setup() \n");
  GradientComputer::setup();
      
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  m_primState.resize(nbEqs);
  
  CFreal dim = static_cast<CFreal>(PhysicalModelStack::getActive()->getDim());
  CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  m_refVol = std::pow(refLength,dim);
  m_refArea = std::pow(refLength,dim-1);
  
  // Connectivity of cells to cell faces
  m_cellFaces =
    Framework::MeshDataStack::getActive()->getConnectivity("cellFaces");

  // To find which TRS the face belongs to
  m_mapGeoToTrs = 
    Framework::MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

  // To build faces
  m_geoBuilder = getMethodData().getFaceTrsGeoBuilder(); 
  Common::SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = m_geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets((*m_globalSockets.getSocketSink<Framework::State*>("states")),
                                (*m_sockets.getSocketSink<Framework::State*>("gstates")), 
                                (*m_globalSockets.getSocketSink<Framework::Node*>("nodes")));
                                
  m_avState.resize(nbEqs);
  m_normal.resize(static_cast<CFuint>(dim));
  
}

//////////////////////////////////////////////////////////////////////////////

CFreal GradientComputerFVMCC::getVolumeAdim(const CFuint& cellID)
{
  DataHandle<CFreal> volumes =  m_sockets.getSocketSink<CFreal>("volumes")->getDataHandle();
  return volumes[cellID]*m_refVol;
}

//////////////////////////////////////////////////////////////////////////////

CFreal GradientComputerFVMCC::getVolume(const CFuint& cellID)
{
  return getVolumeAdim(cellID)*m_refVol;
}

//////////////////////////////////////////////////////////////////////////////

void GradientComputerFVMCC::compute(std::vector<RealVector>& gradients, const CFuint& cellID) 
{
  // Sink sockets required for calculation
  DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();
  DataHandle<CFreal> volumes =  m_sockets.getSocketSink<CFreal>("volumes")->getDataHandle();
  DataHandle<CFreal> normals =  m_sockets.getSocketSink<CFreal>("normals")->getDataHandle();
  DataHandle<CFint> isOutward = m_sockets.getSocketSink<CFint>("isOutward")->getDataHandle();
  DataHandle<RealVector> nodalStates = m_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  
  // reset gradients
  for (CFuint iEq=0; iEq<nbEqs; iEq++) {          
    for (CFuint iDim=0; iDim<dim; ++iDim) {
      gradients[iEq][iDim] = 0.0;
    }
  }
  
  // get the number of faces of element "currStateID"
  const CFuint nbNeighborFaces = m_cellFaces->nbCols(cellID);
  // Loop over all faces of currState
  for (CFuint iFace = 0; iFace < nbNeighborFaces; ++iFace) {
    const CFuint faceID = (*m_cellFaces)(cellID, iFace);
    const CFuint faceIdx = m_mapGeoToTrs->getIdxInTrs(faceID);

    // Get the TRS containing the face
    Common::SafePtr<Framework::TopologicalRegionSet> trs = m_mapGeoToTrs->getTrs(faceID);

    // Construct a face GE
    FaceTrsGeoBuilder::GeoData& geoData = m_geoBuilder->getDataGE();
    geoData.trs = trs;
    geoData.isBFace = m_mapGeoToTrs->isBGeo(faceID);
    geoData.idx = faceIdx;
    m_currFace = m_geoBuilder->buildGE();

    // this is redundant but allows to adapt all the subclasses
    const std::vector<Framework::Node*>* nodesInFace = m_currFace->getNodes();
    const CFuint nbNodesInFace = nodesInFace->size();

   // compute the average state of this face
    m_avState = 0.0;
    for (CFuint i = 0; i < nbNodesInFace; ++i) {
      m_avState += nodalStates[(*nodesInFace)[i]->getLocalID()];
    }
    m_avState /= nbNodesInFace;

    // compute the normal to this face
    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
    for (CFuint iDim=0; iDim<dim; ++iDim)
      m_normal[iDim] = normals[startID+iDim]*m_refArea;

    if (static_cast<CFuint>( isOutward[faceID]) != cellID) 
      m_normal *= -1.;

    // Calculate DIMENSIONAL primitive state
    m_primState = getMethodData().transformToPrimDim(m_avState);

    for (CFuint iEq=0; iEq<nbEqs; iEq++) {          
      for (CFuint iDim=0; iDim<dim; ++iDim) {
        gradients[iEq][iDim] += m_primState[iEq] * m_normal[iDim];
      }
    }

    m_geoBuilder->releaseGE();
  }

  for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
    gradients[iEq] /= (volumes[cellID]*m_refVol);
  }
}


//////////////////////////////////////////////////////////////////////////////

	  } // namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
