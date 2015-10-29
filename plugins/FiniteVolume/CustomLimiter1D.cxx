#include "CustomLimiter1D.hh"
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

MethodStrategyProvider<CustomLimiter1D,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
customLimiter1DProvider("Custom1D");

//////////////////////////////////////////////////////////////////////////////

CustomLimiter1D::CustomLimiter1D(const std::string& name) :
  CustomLimiter(name),
  socket_stencil("stencil"),
  socket_uX("uX")
{
}

//////////////////////////////////////////////////////////////////////////////

CustomLimiter1D::~CustomLimiter1D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter1D::limit(const vector<vector<Node*> >& coord,
			    GeometricEntity* const cell,
			    CFreal* limiterValue)
{
  const GeomEntList *const faces = cell->getNeighborGeos();
  const CFuint nbFaces = faces->size();
  const State *const state = cell->getState(0);
  const CFuint stateID = state->getLocalID();
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();
  DataHandle< vector<State*> > stencil = socket_stencil.getDataHandle();
  
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
      const CFuint nbPoints = 1; // lost generality here
      
      for (CFuint ip = 0; ip < nbPoints; ++ip) {
	computeDeltaMin((*coord[iFace][ip]), *state, iVar);
	
        if (_deltaMin > 0.0) {
	  _gradRatio[0] = deltaPlusMax/_deltaMin; 
	  psi = _functionParser.Eval(&_gradRatio[0]);
	}
	if (_deltaMin < 0.0) {
	  _gradRatio[0] = deltaPlusMin/_deltaMin;
	  psi = _functionParser.Eval(&_gradRatio[0]);
	}
        psimin = min(psi, psimin);
      }
    }
    limiterValue[iVar] = psimin;
    
    if (std::abs(psimin) > 2. || psimin < 0.0) {
      CFout << "CustomLimiter1D::limit() => wrong value for " << _function << " = " << psi << "\n";
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter1D::computeDeltaMin(const Node& coord, const State& state, CFuint iVar)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  const CFuint stateID = state.getLocalID();
  const RealVector& stateCoord = state.getCoordinates();
  _deltaMin = uX(stateID,iVar,state.size())*(coord[XX] - stateCoord[XX]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

