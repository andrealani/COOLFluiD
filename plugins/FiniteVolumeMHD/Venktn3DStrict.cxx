#include "Venktn3DStrict.hh"
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

MethodStrategyProvider<Venktn3DStrict,CellCenterFVMData,Limiter<CellCenterFVMData>,FiniteVolumeModule>
venktn3DStrictProvider("Venktn3DStrict");


void Venktn3DStrict::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
    ("strictCoeff","Fix for smooth flow region.");
  options.addConfigOption< bool >
    ("psiMinEqual1","impose psimin = 1 for certain variables (Haopeng Wang).");
  options.addConfigOption< vector<CFuint> >
	  ("NoLimiterID", "Used to cancel limiter.");
}
//////////////////////////////////////////////////////////////////////////////

Venktn3DStrict::Venktn3DStrict(const std::string& name) :
  Venktn2D(name),
  socket_uZ("uZ")
{
  addConfigOptionsTo(this);
  
  _strictCoeff = 1.0;
  setParameter("strictCoeff",&_strictCoeff);

  _psiMinEqual1 = false;
  setParameter("psiMinEqual1",&_psiMinEqual1);

  _NoLimiterID = vector<CFuint>();
  setParameter("NoLimiterID", &_NoLimiterID);
}

//////////////////////////////////////////////////////////////////////////////

Venktn3DStrict::~Venktn3DStrict()
{
}

//////////////////////////////////////////////////////////////////////////////

void Venktn3DStrict::limit(const vector<vector<Node*> >& coord,
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
 
  CFreal  weight = 0.0;
 
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
    
    //const CFreal deltaPlusMax = max0 - (*state)[iVar];
    //const CFreal deltaPlusMin = min0 - (*state)[iVar];
    CFreal deltaPlusMax = max0 - (*state)[iVar];
    CFreal deltaPlusMin = min0 - (*state)[iVar];
    
    if(iVar==0)
     weight = _strictCoeff*pow(std::abs((min0+1e-10)/(max0+1e-10)),1.0)+(1.0-_strictCoeff);
      
    //weight = 1.0;
    //if((iVar==1 || iVar==2 || iVar==3)&&deltaPlusMax*deltaPlusMin<0.0){ 
      //weight = _strictCoeff*pow(std::abs((min0+1e-20)/(max0+1e-20)),1.0)+(1.0-_strictCoeff);
    //  weight = (std::abs(deltaPlusMin)+1e-20)/(std::abs(deltaPlusMax)+1e-20);
   //   if(weight>1.0) weight=1.0/weight;
   //   weight = _strictCoeff*pow(weight,1.0)+(1.0-_strictCoeff);
   //}
    //if((iVar==1 || iVar==2 || iVar==3)){
    //deltaPlusMax=deltaPlusMax*0.5;
    //deltaPlusMin=deltaPlusMin*0.5;
    //} 
    CFreal psi = 1.0;
    CFreal psimin = 1.1;
    if(_isMFMHD){ // IMPLEMENTATION FROM VENKATAKRISHNAN PAPER 
      const CFreal epsilon2 = pow(_coeffEps*avgDistance/_length, 3.0); //_magnitudeValues[iVar]*_magnitudeValues[iVar]*
      
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	// number of quadrature points associated with this face
	const CFuint nbPoints = 1;
	
	for (CFuint ip = 0; ip < nbPoints; ++ip) {
	  Venktn3DStrict::computeDeltaMin((*coord[iFace][ip]), *state, iVar);
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
          if(iVar<4 || iVar==7)
	    psimin = min(psi*weight, psimin);
          else
	   psimin = (!_psiMinEqual1) ? min(psi, psimin) : 1.0;
	}
      }
    
     /* for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const GeometricEntity *const neighborFace = (*faces)[iFace];
	State *const neighborState = (state == neighborFace->getState(0)) ?
	  neighborFace->getState(1) : neighborFace->getState(0);
	
	      Venktn3DStrict::computeDeltaMin(neighborState->getCoordinates(), *state, iVar);
              const CFreal dU = 0.8*(*neighborState)[iVar] - (*state)[iVar];
              if(_deltaMin * dU < 0.0){
                psi = 0.0;
              }
	      else if (_deltaMin > 0.0) {
	        const CFreal dMinStar = MathTools::MathFunctions::sign(_deltaMin)*(std::abs(_deltaMin)/_magnitudeValues[iVar] + 1e-24);
	        const CFreal dPlus = dU/_magnitudeValues[iVar];
	        const CFreal dPlus2   = dPlus*dPlus;
	        const CFreal dPlusMin = dPlus*dMinStar;
	        const CFreal Num = (dPlus2 + 0.001*epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
	        const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + 0.001*epsilon2);
	        psi = 1./dMinStar*(Num/Den);
	       // psi = min(1.0, dU/_deltaMin);
	      }
	      else if (_deltaMin < 0.0) {
	        const CFreal dMinStar = MathTools::MathFunctions::sign(_deltaMin)*(std::abs(_deltaMin)/_magnitudeValues[iVar] + 1e-24);
	        const CFreal dPlus = dU/_magnitudeValues[iVar];
	        const CFreal dPlus2   = dPlus*dPlus;
	        const CFreal dPlusMin = dPlus*dMinStar;
	        const CFreal Num = (dPlus2 + 0.001*epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
	        const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + 0.001*epsilon2);
	        psi = 1./dMinStar*(Num/Den);
	       // psi = min(1.0, dU/_deltaMin);
	       }
	  psimin2 = min(psi, psimin2);
              }
      psimin = weight * psimin + (1.0 - weight) * psimin2; */ 
    }
    limiterValue[iVar] = psimin;
    
    const CFreal maxAllowableLimiterFunctionValue = 1.094;
    if (limiterValue[iVar] > maxAllowableLimiterFunctionValue) {
      CFout << "wrong limiterValue = " << limiterValue[iVar] << "\n";
    }
  }
  //>> Mark 2025.07.30 by Hp to cancel limiter for several reconstructed variables
  if (_NoLimiterID.size() > 0){
	  for (CFuint iVar = 0; iVar < _NoLimiterID.size(); ++iVar){
		  limiterValue[_NoLimiterID[iVar]] = 1.0;
	  }
  }
  //<< Mark 2025.07.30 by Hp to cancel limiter for several reconstructed variables
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

