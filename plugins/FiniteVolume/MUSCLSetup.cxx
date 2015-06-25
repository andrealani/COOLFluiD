#include "FiniteVolume/FiniteVolume.hh"
#include "MUSCLSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MUSCLSetup, 
		      CellCenterFVMData, 
		      FiniteVolumeModule> 
upwindBiasedMUSCLSetupProvider("MUSCLSetup");

//////////////////////////////////////////////////////////////////////////////

MUSCLSetup::MUSCLSetup(const std::string& name) :
  StdSetup(name),
  socket_stencil("stencil"),
  socket_uX("uX"),
  socket_uY("uY"),
  socket_uZ("uZ")
{
}

//////////////////////////////////////////////////////////////////////////////

MUSCLSetup::~MUSCLSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLSetup::configure ( Config::ConfigArgs& args )
{
  StdSetup::configure(args);
  
  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);
  
  // set the number of overlap layers to two for parallel runs
  if (PE::GetPE().IsParallel()) {
    meshData->setNbOverlapLayers(3);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MUSCLSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
    StdSetup::providesSockets();
  result.push_back(&socket_stencil);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLSetup::execute()
{
  CFAUTOTRACE;
  
  StdSetup::execute();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();
  cf_assert(nbStates > 0);
  
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  cf_assert(nbFaces > 0);
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  stencil.resize(nbFaces);
    
  GeometricEntityPool<FaceTrsGeoBuilder> faceBuilder;
  faceBuilder.setup();
  faceBuilder.getGeoBuilder()->setDataSockets
    (socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = faceBuilder.getDataGE();
  
  GeometricEntityPool<CellTrsGeoBuilder> cellBuilder;
  cellBuilder.setup();
  cellBuilder.getGeoBuilder()->setDataSockets
    (socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder.getDataGE();
  
  SafePtr<TopologicalRegionSet> cells = 
    MeshDataStack::getActive()->getTrs("InnerCells");
  cellData.trs = cells;
  
  // set the list of faces
  vector<Common::SafePtr<TopologicalRegionSet> > trs = 
    MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();
  
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    faceData.trs = currTrs;
    
    // the faces on the boundary of the partition don't have to
    // be processed (their fluxes could give NaN)
    if ((!currTrs->hasTag("partition")) && currTrs->hasTag("face")) {
      if (currTrs->hasTag("writable")) {
	faceData.isBFace = true;
      }
      else {
        faceData.isBFace = false;
      }
      
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        CFLogDebugMed( "iFace = " << iFace << "\n");
	
        // build the GeometricEntity
        faceData.idx = iFace;
        GeometricEntity *const face = faceBuilder.buildGE();
	const CFuint faceID = face->getID();
	stencil[faceID].reserve((faceData.isBFace == false) ? 4 : 3);
	stencil[faceID].push_back(face->getState(LEFT));
	stencil[faceID].push_back(face->getState(RIGHT));
	
	// left cell
	cellData.idx = face->getState(LEFT)->getLocalID();
	
	CFuint oppositeIFace = 0;
	GeometricEntity* cell = cellBuilder.buildGE();
	const vector<GeometricEntity*>& cellFacesL = *cell->getNeighborGeos();
	const CFuint nbCellFacesL =  cellFacesL.size();
	for (CFuint f = 0; f < nbCellFacesL; ++f) {
	  if (cellFacesL[f]->getID() == faceID) {
	    oppositeIFace = getMethodData().getOppositeIFace
	      (f, PhysicalModelStack::getActive()->getDim(), cell->nbNodes());
	    
	    State *const leftState  = cellFacesL[oppositeIFace]->getState(0);
	    State *const rightState = cellFacesL[oppositeIFace]->getState(1);
	    stencil[faceID].push_back((cell->getState(0) == leftState) ? 
				      rightState : leftState);
	    break;
	  }
	}
	cellBuilder.releaseGE();
	
	// if the current face is not a boundary face, look for the 
	// distance-1 neighbor of the face
	if (faceData.isBFace == false) {
	  cellData.idx = face->getState(RIGHT)->getLocalID();
	  cell = cellBuilder.buildGE();
	  	  
	  const vector<GeometricEntity*>& cellFacesR = *cell->getNeighborGeos();
	  const CFuint nbCellFacesR =  cellFacesR.size();
	  for (CFuint f = 0; f < nbCellFacesR; ++f) {
	    if (cellFacesR[f]->getID() == faceID) {
	      oppositeIFace = getMethodData().getOppositeIFace
		(f, PhysicalModelStack::getActive()->getDim(), cell->nbNodes());
	      
	      State *const leftState  = cellFacesR[oppositeIFace]->getState(0);
	      State *const rightState = cellFacesR[oppositeIFace]->getState(1);
	      stencil[faceID].push_back((cell->getState(0) == leftState) ? 
					rightState : leftState);
	      break;
	    }
	  }
	  cellBuilder.releaseGE();
	}
	
	faceBuilder.releaseGE();
	
	// // sanity check (could become optional) 
// 	set<State *> setS;
// 	const CFuint sizeS = stencil[faceID].size();
// 	for (CFuint j = 0; j < sizeS; ++j) {
// 	  setS.insert(stencil[faceID][j]);  
// 	}
// 	if (setS.size() != sizeS)  {
// 	  cout << sizeS << " => ERROR STENCIL on face " << faceID << " ";
// 	  for (CFuint j = 0; j < sizeS; ++j) {
// 	    cout << stencil[faceID][j] << " "; 
// 	  }
// 	  cout << currTrs->getName() << ", iTRS = " << iTRS << ", iFace = " << iFace << endl;
// 	  abort();
// 	}
      } 
    }
  }
  
  // resize the limiter storage
  // by default would have size == 0
  // in this case the limiter is facewise and not cellwise as in the LeastSquareP1Setup
  const CFuint sizeLimiter = 2*nbFaces*PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> limiter = socket_limiter.getDataHandle();
  limiter.resize(sizeLimiter);
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  uX.resize(nbStates*nbEqs);
  
  if (PhysicalModelStack::getActive()->getDim() >= DIM_2D) {
    DataHandle<CFreal> uY = socket_uY.getDataHandle();
    uY.resize(nbStates*nbEqs);
  }
  
  if (PhysicalModelStack::getActive()->getDim() == DIM_3D) {
    DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
    uZ.resize(nbStates*nbEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
