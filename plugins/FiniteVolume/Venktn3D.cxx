#include "Venktn3D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "MathTools/MathFunctions.hh"
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

MethodStrategyProvider<Venktn3D,CellCenterFVMData,Limiter<CellCenterFVMData>,FiniteVolumeModule>
venktn3DProvider("Venktn3D");

//////////////////////////////////////////////////////////////////////////////

Venktn3D::Venktn3D(const std::string& name) :
  Venktn2D(name),
  socket_uZ("uZ")
{
}

//////////////////////////////////////////////////////////////////////////////

Venktn3D::~Venktn3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Venktn3D::limit(const vector<vector<Node*> >& coord,
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
  
  for(CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    CFreal min0 = (*state)[iVar];
    CFreal max0 = (*state)[iVar];
    CFreal avgDistance = 0.0;
    
    if (this->m_useFullStencil) {
      CFuint nbNeighbors = 0;
      if (!this->m_useNodalExtrapolationStencil) {
	const vector<State*>& s = stencil[stateID];
	nbNeighbors = s.size();
	for (CFuint is = 0; is < nbNeighbors; ++is) {
	  State *const currState = s[is];
	  if(state != currState) {
	    min0 = min(min0,(*currState)[iVar]);
	    max0 = max(max0,(*currState)[iVar]);
	    avgDistance += MathFunctions::getDistance
	      (state->getCoordinates(), currState->getCoordinates());
	  }
	}
      }
      else {
	// old method kept for backward compatibility: it may escludes some ghost states at corners
	const CFuint nbNodesInCell = cell->nbNodes();
	for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
	  const CFuint nodeID = cell->getNode(iNode)->getLocalID();
	  const vector<State*>& s = getMethodData().
	    getNodalStatesExtrapolator()->getNodalStateNeighbors(nodeID);
	  for (CFuint is = 0; is < s.size(); ++is) {
	    State *const currState = s[is];
	    if(state != currState) {
	      min0 = min(min0,(*currState)[iVar]);
	      max0 = max(max0,(*currState)[iVar]);
	      avgDistance += MathFunctions::getDistance
		(state->getCoordinates(), currState->getCoordinates());
	      nbNeighbors++;
	    }
	  }
	}
      }
      
      // only uniquely counted neighbors should be considered ...
      avgDistance /= nbNeighbors;
    }
    else {
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const GeometricEntity *const neighborFace = (*faces)[iFace];
	State *const neighborState = (state == neighborFace->getState(0)) ?
	  neighborFace->getState(1) : neighborFace->getState(0);
	
	min0 = min(min0,(*neighborState)[iVar]);
	max0 = max(max0,(*neighborState)[iVar]);
	
	// the average distance between the centroids of the cell and
	// each of its face neighbors is calculated
	avgDistance += MathFunctions::getDistance
	  (state->getCoordinates(), neighborState->getCoordinates());
      }
      
      avgDistance /= nbFaces;
    }
    
    const CFreal deltaPlusMax = max0 - (*state)[iVar];
    const CFreal deltaPlusMin = min0 - (*state)[iVar];
    CFreal psi = 1.0;
    CFreal psimin = 1.1;
    
    if(_isMFMHD){ // IMPLEMENTATION FROM VENKATAKRISHNAN PAPER 
      const CFreal epsilon2 = pow(_coeffEps*avgDistance/_length, 3.0); //_magnitudeValues[iVar]*_magnitudeValues[iVar]*
      
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	// number of quadrature points associated with this face
	const CFuint nbPoints = 1;
	
	for (CFuint ip = 0; ip < nbPoints; ++ip) {
	  Venktn3D::computeDeltaMin((*coord[iFace][ip]), *state, iVar);
	  
	  if (_deltaMin > 0.0) {
	    const CFreal dMinStar = MathTools::MathFunctions::sign(_deltaMin)*(std::abs(_deltaMin)/_magnitudeValues[iVar] + 1e-24);
	    const CFreal dPlus = deltaPlusMax/_magnitudeValues[iVar];
	    const CFreal dPlus2   = dPlus*dPlus;
	    const CFreal dPlusMin = dPlus*dMinStar;
	    const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
	    const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
	    psi = 1./dMinStar*(Num/Den);
	  }
	  if (_deltaMin < 0.0) {
	    const CFreal dMinStar = MathTools::MathFunctions::sign(_deltaMin)*(std::abs(_deltaMin)/_magnitudeValues[iVar] + 1e-24);
	    const CFreal dPlus = deltaPlusMin/_magnitudeValues[iVar];
	    const CFreal dPlus2   = dPlus*dPlus;
	    const CFreal dPlusMin = dPlus*dMinStar;
	    const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
	    const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
	    psi = 1./dMinStar*(Num/Den);
	  }
	  psimin = min(psi, psimin);
	}
      }
    }
    else{ //OLD IMPLEMENTATION
      const CFreal epsilon = _coeffEps*_magnitudeValues[iVar]*_magnitudeValues[iVar]*
	pow(avgDistance/_length, 3.0);
      
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	// number of quadrature points associated with this face
	const CFuint nbPoints = 1;
	
	for (CFuint ip = 0; ip < nbPoints; ++ip) {
	  const CFreal xq = (*coord[iFace][ip])[0];
	  const CFreal yq = (*coord[iFace][ip])[1];
	  const CFreal zq = (*coord[iFace][ip])[2];
	  const CFreal deltaMin = uX(stateID,iVar,nbEquations)*
	    (xq - stateCoord[0]) + uY(stateID,iVar,nbEquations)*
	    (yq - stateCoord[1]) + uZ(stateID,iVar,nbEquations)*
	    (zq - stateCoord[2]);
	  
	  if (deltaMin > 0.0) {
	    const CFreal dPlusMax2   = deltaPlusMax*deltaPlusMax;
	    const CFreal dPlusMaxMin = deltaPlusMax*deltaMin;
	    // epsilon is a scaling factor which is assigned arbitrarily
	    psi = (dPlusMax2 + 2.0*dPlusMaxMin + epsilon)/
	      (dPlusMax2 + dPlusMaxMin + 2.0*deltaMin*deltaMin + epsilon);
	  }
	  if (deltaMin < 0.0) {
	    const CFreal y = deltaPlusMin/deltaMin;
	    psi = (y*y + 2*y)/(y*y + y + 2.);
	  }
	  psimin = min(psi, psimin);
	}
      }
    }
    limiterValue[iVar] = psimin;
    
    const CFreal maxAllowableLimiterFunctionValue = 1.094;
    if (limiterValue[iVar] > maxAllowableLimiterFunctionValue) {
      CFout << "wrong limiterValue = " << limiterValue[iVar] << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

