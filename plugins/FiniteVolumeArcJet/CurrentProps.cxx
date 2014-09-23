#include "Common/PE.hh"
#include "MathTools/MathConsts.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Common/BadValueException.hh"
#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeArcJet/CurrentProps.hh"


/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CurrentProps, DataProcessingData, FiniteVolumeArcJetModule>
CurrentPropsProvider("CurrentProps");

//////////////////////////////////////////////////////////////////////////////

void CurrentProps::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

CurrentProps::CurrentProps(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_sigma("sigma"),
  socket_Jx("Jx"),
  socket_Jy("Jy"),
  socket_Jz("Jz"),
  m_gradPhi(),
  m_normal()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

CurrentProps::~CurrentProps()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
CurrentProps::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  
  result.push_back(&socket_sigma);
  result.push_back(&socket_Jx);
  result.push_back(&socket_Jy);
  result.push_back(&socket_Jz);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CurrentProps::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);  
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CurrentProps::setup()
{
  CFAUTOTRACE;

  CFout << "Setup: Computing electric current properties.\n";
  
  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();
  const CFuint DIM = 3; //AAL: to generalize in the future for 2D

  DataHandle<CFreal> sigma = socket_sigma.getDataHandle();
  sigma.resize(nbCells);
  sigma = 0.0;
  
  DataHandle<CFreal> Jx = socket_Jx.getDataHandle();
  Jx.resize(nbCells);
  Jx = 0.0;
  
  DataHandle<CFreal> Jy = socket_Jy.getDataHandle();
  Jy.resize(nbCells);
  Jy = 0.0;
  
  DataHandle<CFreal> Jz = socket_Jz.getDataHandle();
  Jz.resize(nbCells);
  Jz = 0.0;

  m_geoBuilder.setup();


  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(m_library.isNotNull());
  
  m_gradPhi.resize(DIM, 0.);
  m_normal.resize(DIM, 0.);
}

//////////////////////////////////////////////////////////////////////////////

void CurrentProps::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> sigma = socket_sigma.getDataHandle();
  DataHandle<CFreal> Jx = socket_Jx.getDataHandle();
  DataHandle<CFreal> Jy = socket_Jy.getDataHandle();
  DataHandle<CFreal> Jz = socket_Jz.getDataHandle();
  
  const CFuint DIM = 3; //AAL: to generalize to 2D in the future
  
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  //Comment to know the name of the TRs
  //for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs) {
    //SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
    //cout <<"List of Trs" << trs->getName() <<"\n";
  //}
  

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint nbCells = cells->getLocalNbGeoEnts();

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    sigma[iCell] = 0;
    m_gradPhi = 0.;
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID();
    
    const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();    
    const CFuint phiID = totalNbEqs-1;
    const CFuint TID = phiID-1;
    const CFuint pID = 0;
    
    State *currState = currCell->getState(0);
   
    CFreal* tVec = CFNULL;
    CFreal Tdim = (*currState)[TID];
    CFreal pdim = (*currState)[pID];
    const CFreal sigmaLocal = m_library->sigma(Tdim, pdim, tVec);
    
    sigma[iCell] = sigmaLocal;
    
    const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
    const CFuint nbFaces = faces.size(); 
    const CFreal ovVolume = 1./volumes[elemID];
//     cout << "faces"<< &faces <<"\n";
//     cout << "nbFaces"<< nbFaces <<"\n"; 
    for (CFuint i = 0; i < nbFaces; ++i) {
      const GeometricEntity *const face = currCell->getNeighborGeo(i);
      const CFuint faceID = face->getID();
      const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
      const CFuint nbFaceNodes = face->nbNodes();
      const CFreal ovNbFaceNodes = 1./(CFreal)nbFaceNodes;
      //cout << "Calculating gradPhi 2\n";
      
      // store the outward face normal (scaled with the corresponding face area)
      // fill in only the components that are available (==dim)
      // in 2D, m_normal[ZZ] = 0
      const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
      for (CFuint d = 0; d < DIM; ++d) {
	m_normal[d] = normals[startID+d]*factor;
      }
      //cout << "Calculating gradPhi 3\n";
//    Comments, to be removed in the future in order to include the effect of UxB   
//       // compute contribution of this face to the Green Gauss integral:
//       // \bar{\vec{B}}_f \vec{n}_f = \sum_{i \in f} B_i/N_f \vec{n}_f
//       m_uf = 0.;
//       m_xf = 0.;
//       CFreal avTdim = 0.;
//       CFreal avpdim = 0.;
//       CFreal* tVec = CFNULL;
//       cout << "Before the faces\n";   
      for (CFuint n = 0; n < nbFaceNodes; ++n) {
	const CFuint nodeID = face->getNode(n)->getLocalID();
	const RealVector& nodalState = nstates[nodeID];
// 	// consider all 3D components here even in 2D
	for (CFuint d = 0; d < DIM; ++d) {
	  const CFreal nd = m_normal[d];
	  m_gradPhi[d] += nd*nodalState[phiID]*ovNbFaceNodes;
// 	  cout << "nodalState[phiID]" << nodalState[phiID] << "\n";
// 	  m_uf[d]      += nodalState[1+d]; // change here for NEQ
	}
	//cout << "Calculating gradPhi 4\n";
// 	
// 	m_xf   += *face->getNode(n);
// 	avTdim += nodalState[TID]; 
// 	avpdim += nodalState[pID];
      }
//       
//       m_uf   *= ovNbFaceNodes;
//       m_xf   *= ovNbFaceNodes;
//       avTdim *= ovNbFaceNodes; 
//       avpdim *= ovNbFaceNodes;
//       
//       // AL: compute m_Bf in function of the midFace point
//       //     Rough assumptions here ...
//       computeB(m_xf[XX] - m_electrodeX, m_Bf);
//       const CFreal sigma = m_library->sigma(avTdim, avpdim, tVec);
// //       cout <<"ArcJetPhiST::computeSource => sigma =" << sigma <<"\n";
//       MathFunctions::crossProd(m_uf, m_Bf, m_UxB);
//       sigmaUxB += sigma*MathFunctions::innerProd(m_UxB, m_normal);
    }
//     cout << "gradPhi" << m_gradPhi << "\n";
   
    m_gradPhi *= ovVolume;
    
    Jx[iCell] = -sigmaLocal*m_gradPhi[0];
    Jy[iCell] = -sigmaLocal*m_gradPhi[1];
    Jz[iCell] = -sigmaLocal*m_gradPhi[2];
    
    // Set the boundary at 0
    
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
      State *const neighborState = (currCell->getState(0) == face->getState(0)) ? face->getState(1) : face->getState(0);
      if (neighborState->isGhost()){
	Jx[iCell] = 0;
	Jy[iCell] = 0;
	Jz[iCell] = 0;
      }
    }
    
    m_geoBuilder.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void CurrentProps::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

