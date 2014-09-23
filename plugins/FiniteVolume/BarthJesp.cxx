#include "BarthJesp.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Common/CFLog.hh"
#include "FiniteVolume/FiniteVolume.hh"

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

MethodStrategyProvider<BarthJesp,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
barthJespProvider("BarthJesp");
      
MethodStrategyProvider<BarthJesp,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
barthJesp2DProvider("BarthJesp2D");

MethodStrategyProvider<BarthJesp,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
barthJesp3DProvider("BarthJesp3D");
      
//////////////////////////////////////////////////////////////////////////////
      
BarthJesp::BarthJesp(const std::string& name) :
  Limiter<CellCenterFVMData>(name),
  socket_stencil("stencil"),
  socket_uX("uX"),
  socket_uY("uY"),
  socket_uZ("uZ")
{
  // the following default value is to keep compatibility with old testcases
  // user can still choose to override this value in the CFcase file
  if (getName() == "BarthJesp3D") {
    m_useFullStencil = false;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

BarthJesp::~BarthJesp()
{
}

//////////////////////////////////////////////////////////////////////////////

void BarthJesp::limit(const vector<vector<Node*> >& coord,
		      GeometricEntity* const cell,
		      CFreal* limiterValue)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  DataHandle< vector<State*> > stencil = socket_stencil.getDataHandle();
  
  const GeomEntList *const faces = cell->getNeighborGeos();
  const CFuint nbFaces = faces->size();
  const State *const state = cell->getState(0);
  const CFuint stateID = state->getLocalID();
  const RealVector& stateCoord = state->getCoordinates();
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  for(CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    CFreal min0 = (*state)[iVar];
    CFreal max0 = (*state)[iVar];
    
    if (this->m_useFullStencil) { 
      if (!this->m_useNodalExtrapolationStencil) {
	const vector<State*>& s = stencil[stateID];
	for (CFuint is = 0; is < s.size(); ++is) {
	  State *const currState = s[is];
	  if(state != currState) {
	    min0 = min(min0,(*currState)[iVar]);
	    max0 = max(max0,(*currState)[iVar]);
	  }
	}
      }
      else {
	// old method kept for backward compatibility: it may escludes some ghost states at corners
	const CFuint nbNodesInCell = cell->nbNodes();
	for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
	  const CFuint nodeID = cell->getNode(iNode)->getLocalID();
	  const vector<State*>& s = getMethodData().getNodalStatesExtrapolator()->getNodalStateNeighbors(nodeID);
	  for (CFuint is = 0; is < s.size(); ++is) {
	    State *const currState = s[is];
	    if(state != currState) {
	      min0 = min(min0,(*currState)[iVar]);
	      max0 = max(max0,(*currState)[iVar]);
	    }
	  }
	}
      }
    }
    else {
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const GeometricEntity *const neighborFace = (*faces)[iFace];
	State *const neighborState = (state == neighborFace->getState(0)) ?
	  neighborFace->getState(1) : neighborFace->getState(0);
	min0 = min(min0,(*neighborState)[iVar]);
	max0 = max(max0,(*neighborState)[iVar]);
      }
    }
    
    const CFreal deltaPlusMax = max0 - (*state)[iVar];
    const CFreal deltaPlusMin = min0 - (*state)[iVar];
    CFreal psi = 1.0;
    CFreal psimin = 1.1;
    
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      
      // number of quadrature points associated with this face
      const CFuint nbPoints = 1; // here somebody should give this info
      for (CFuint ip = 0; ip < nbPoints; ++ip) {
        const CFreal xq = (*coord[iFace][ip])[XX];
        const CFreal yq = (*coord[iFace][ip])[YY];
	CFreal deltaMin = uX(stateID,iVar,nbEquations)*
          (xq - stateCoord[XX]) + uY(stateID,iVar,nbEquations)*
          (yq - stateCoord[YY]);
        if (dim == DIM_3D) {
	  const CFreal zq = (*coord[iFace][ip])[ZZ];
	  deltaMin += uZ(stateID,iVar,nbEquations)*(zq - stateCoord[ZZ]);
	}
	
	if (deltaMin > 0.0) {
          psi = min((CFreal)1.0, deltaPlusMax/deltaMin);
        }
        if (deltaMin < 0.0) {
          psi = min((CFreal)1.0, deltaPlusMin/deltaMin);
        }
        psimin = min(psi, psimin);
      }
    }
    
    limiterValue[iVar] = psimin;
    
    if (limiterValue[iVar] > 1.0) {
      CFout << "wrong limiterValue = " << limiterValue[iVar] << "\n";
    }
  }
} 
      
//////////////////////////////////////////////////////////////////////////////
  
std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > BarthJesp::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
    Framework::Limiter<CellCenterFVMData>::needsSockets();
  result.push_back(&socket_stencil);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);  
  result.push_back(&socket_uZ);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

