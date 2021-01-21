#include "LeastSquareP1PolyRec2D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

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

MethodStrategyProvider<LeastSquareP1PolyRec2D, CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeModule> 
leastSquareP1PolyRec2DProvider("LinearLS2D");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2D::LeastSquareP1PolyRec2D(const std::string& name) :
  FVMCC_PolyRec(name),
  socket_stencil("stencil"),
  socket_weights("weights"),
  socket_uX("uX"),
  socket_uY("uY"),
  _l11(),
  _l12(),
  _l22(),
  _lf1(),
  _lf2()
{
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2D::~LeastSquareP1PolyRec2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2D::configure ( Config::ConfigArgs& args )
{
  FVMCC_PolyRec::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LeastSquareP1PolyRec2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    FVMCC_PolyRec::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_stencil);
  result.push_back(&socket_weights);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
      
void LeastSquareP1PolyRec2D::computeGradients()
{
  CFLog(VERBOSE, "LeastSquareP1PolyRec2D::computeGradients() => START\n");
  
  prepareReconstruction();
 
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> weights = socket_weights.getDataHandle();
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();

  for(CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    _lf1 = 0.0;
    _lf2 = 0.0;
    CFuint iEdge = 0;
    
    for(CFuint iState = 0; iState < nbStates; ++iState) {
      assert(iState == states[iState]->getLocalID());
      const State* const first = states[iState];
      const CFuint stencilSize = stencil[iState].size();
      // loop over the neighbor cells belonging to the chosen stencil
      for(CFuint in = 0; in < stencilSize; ++in) {
	const State* const last = stencil[iState][in];
	const CFuint lastID = (!last->isGhost()) ? last->getLocalID() : 
	  numeric_limits<CFuint>::max();
	const CFuint firstID = first->getLocalID();
	cf_assert(firstID != lastID);
	
	if (lastID > firstID) {
	  // consider the next edge
	  const RealVector& nodeFirst = first->getCoordinates();
	  const RealVector& nodeLast = last->getCoordinates();
	  const CFreal weig = weights[iEdge];
	  const CFreal dx = weig*(nodeLast[0] - nodeFirst[0]);
	  const CFreal dy = weig*(nodeLast[1] - nodeFirst[1]);
	  const CFreal du = weig*((*last)[iVar] - (*first)[iVar]);
	  const CFreal dxdu = dx*du;
	  const CFreal dydu = dy*du;
	  
	  _lf1[firstID] += dxdu;
	  _lf2[firstID] += dydu;
	  
	  if (!last->isGhost()) {
	    _lf1[lastID] += dxdu;
	    _lf2[lastID] += dydu;
	  }
	  
	  ++iEdge;
	}
      }
    }
    
    for(CFuint iState = 0; iState < nbStates; ++iState) {
      const CFreal invDet = 1./(_l11[iState]*_l22[iState] - _l12[iState]*_l12[iState]);
      uX(iState,iVar,nbEquations) = (_l22[iState]*_lf1[iState] - _l12[iState]*_lf2[iState])*invDet;
      uY(iState,iVar,nbEquations) = (_l11[iState]*_lf2[iState] - _l12[iState]*_lf1[iState])*invDet;
      
      CFLogDebugMax( "invDet = " << invDet
		     << ", l11 = " << _l11[iState]
		     << ", l12 = " << _l12[iState]
		     << ", l22 = " << _l22[iState]
		     <<", lf1 = " << _lf1[iState]
		     <<", lf2 = " << _lf2[iState]
		     << ", uX =" << uX(iState,iVar,nbEquations)
		     << ", uY =" << uY(iState,iVar,nbEquations)
		     << "\n");
    }
  }
  
  CFLog(VERBOSE, "LeastSquareP1PolyRec2D::computeGradients() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2D::extrapolateImpl(GeometricEntity* const face)
{
  FVMCC_PolyRec::baseExtrapolateImpl(face);
  
  // first compute the limiter
  computeFaceLimiter(face);
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  // please note that you work with references !!!
  // don't forget the "&" !!!
  const vector<Node*>& coord = getCoord();
  vector<State*>& valuesLR  = getExtrapolatedValues();
  const State *const state = face->getState(LEFT);
  const CFuint stateID = state->getLocalID();
  const RealVector& stateCoord = state->getCoordinates();
  const State *const neighState = face->getState(RIGHT);
  const CFuint neighStateIDgrad = (!isBoundaryFace()) ? neighState->getLocalID() : stateID;
  const RealVector& neighStateCoord = neighState->getCoordinates();

  // number of quadrature points associated with this face
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint startID = stateID*nbEqs;
  const CFuint neighStartID = neighStateIDgrad*nbEqs;
  cf_assert(_zeroGradient != CFNULL);
  
  const CFreal xq = (*coord[0])[XX];
  const CFreal yq = (*coord[0])[YY];
  
  if (_vars.size() > 1) {
    _vFunction.evaluate((*coord[0]),_gradientCoeff);
  }
  
  for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    cf_assert(iVar < _zeroGradient->size());
    
    // AL: gradientCoeff prevails on the interactive gradientFactor
    if (_vars.size() > 1) {_gradientFactor = _gradientCoeff[iVar];}
    
    const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ? 
      _gradientFactor*(uX(stateID,iVar,nbEqs)*(xq - stateCoord[0]) + 
		       uY(stateID,iVar,nbEqs)*(yq - stateCoord[1])) : 0.0;
    
    const CFreal gradientCoeffNeighbor = (!(*_zeroGradient)[iVar]) ? 
      _gradientFactor*(uX(neighStateIDgrad,iVar,nbEqs)*(xq - neighStateCoord[0]) + 
		       uY(neighStateIDgrad,iVar,nbEqs)*(yq - neighStateCoord[1])) : 0.0;
    
    // left reconstructed value (the one inside the current cell)
    const CFreal lValue = (_stopLimiting == 0) ? newLimiter[startID + iVar] : 1.;
    (*valuesLR[0])[iVar] = (*state)[iVar] + lValue*gradientCoeffState;
    
    // right reconstructed value (the one inside the neighbor cell)
    const CFreal rValue = (_stopLimiting == 0) ? newLimiter[neighStartID + iVar] : 1.;
    (*valuesLR[1])[iVar] = (*neighState)[iVar] + rValue*gradientCoeffNeighbor;
    
    getBackupValues(0)[iVar] =  (*valuesLR[0])[iVar];
    getBackupValues(1)[iVar] =  (*valuesLR[1])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2D::extrapolateImpl(GeometricEntity* const face,
					     CFuint iVar, CFuint leftOrRight)
{
  CFLog(DEBUG_MED, "LeastSquareP1PolyRec2D::extrapolateImpl() => start\n");
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
  const State *const state = face->getState(leftOrRight);
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFuint stateID = 0;
  if (leftOrRight == LEFT) {
    stateID = state->getLocalID();
  }
  else {
    stateID = (!isBoundaryFace()) ? state->getLocalID() :
      face->getState(LEFT)->getLocalID();
  }
  const RealVector& stateCoord = state->getCoordinates();
  
  // number of quadrature points associated with this face
  const vector<Node*>& coord = getCoord();
  const CFuint startID = stateID*nbEqs;
  
  const CFreal xq = (*coord[0])[XX];
  const CFreal yq = (*coord[0])[YY];
  copyBackupValues();
  
  if (_vars.size() > 1) {
    _vFunction.evaluate(iVar,(*coord[0]), _gradientFactor);
  }
  
  // AL: gradientCoeff prevails on the interactive gradientFactor
  const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ?
    _gradientFactor*(uX(stateID,iVar,nbEqs)*(xq - stateCoord[XX]) + uY(stateID,iVar,nbEqs)*(yq - stateCoord[YY])) : 0.0;
  
  // reconstructed value
  const CFreal lValue = (_stopLimiting == 0) ? newLimiter[startID + iVar] : 1.;
  getValues(leftOrRight)[iVar] = (*state)[iVar] + lValue*gradientCoeffState;
  
  CFLog(DEBUG_MED, "LeastSquareP1PolyRec2D::extrapolateImpl() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2D::updateWeights()
{

  FVMCC_PolyRec::updateWeights();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> weights = socket_weights.getDataHandle();

  const CFuint nbStates = states.size();

  _l11 = 0.0;
  _l12 = 0.0;
  _l22 = 0.0;
  _lf1 = 0.0;
  _lf2 = 0.0;

  // weight coefficients are calculated
 CFuint iEdge = 0;
 for(CFuint iState = 0; iState < nbStates; ++iState) {
   const State* const first = states[iState];
   const CFuint stencilSize = stencil[iState].size();

   // loop over the neighbor cells belonging to the chosen stencil
   for(CFuint in = 0; in < stencilSize; ++in) {
     const State* const last = stencil[iState][in];
     const CFuint lastID = (!last->isGhost()) ? last->getLocalID() : 
       numeric_limits<CFuint>::max();
     const CFuint firstID = first->getLocalID();
     cf_assert(firstID != lastID);
     
     if (lastID > firstID) {
       const RealVector& nodeFirst = first->getCoordinates();
       const RealVector& nodeLast = last->getCoordinates();
       const CFreal deltaR = MathFunctions::getDistance(nodeFirst,nodeLast);
       
       weights[iEdge] = 1.0/deltaR;

       // weights always != 0
       const CFreal dx = weights[iEdge]*(nodeLast[0]
            - nodeFirst[0]);
       const CFreal dy = weights[iEdge]*(nodeLast[1]
            - nodeFirst[1]);

       CFLogDebugMax( "weights = " << weights[iEdge]
       << "dx = " << dx
       << "dy = " << dy << "\n");

       _l11[firstID] += dx*dx;
       _l12[firstID] += dx*dy;
       _l22[firstID] += dy*dy;
       
       if (!last->isGhost()) {
	 _l11[lastID] += dx*dx;
	 _l12[lastID] += dx*dy;
	 _l22[lastID] += dy*dy;
       }
       
       ++iEdge;
     }
   }
 }

}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2D::setup()
{
  FVMCC_PolyRec::setup();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> weights = socket_weights.getDataHandle();

  const CFuint nbStates = states.size();

  _l11.resize(nbStates);
  _l12.resize(nbStates);
  _l22.resize(nbStates);
  _lf1.resize(nbStates);
  _lf2.resize(nbStates);

  _l11 = 0.0;
  _l12 = 0.0;
  _l22 = 0.0;
  _lf1 = 0.0;
  _lf2 = 0.0;

  // weight coefficients are calculated
  CFuint iEdge = 0;
  for(CFuint iState = 0; iState < nbStates; ++iState) {
   const State* const first = states[iState];
   const CFuint stencilSize = stencil[iState].size();

   // loop over the neighbor cells belonging to the chosen stencil
   for(CFuint in = 0; in < stencilSize; ++in) {
     const State* const last = stencil[iState][in];
     const CFuint lastID = (!last->isGhost()) ? last->getLocalID() : 
       numeric_limits<CFuint>::max();
     const CFuint firstID = first->getLocalID();
     cf_assert(firstID != lastID);
     
     if (lastID > firstID) {
       const RealVector& nodeFirst = first->getCoordinates();
       const RealVector& nodeLast = last->getCoordinates();
       const CFreal deltaR = MathFunctions::getDistance(nodeFirst,nodeLast);

       weights[iEdge] = 1.0/deltaR;

       // weights always != 0
       const CFreal dx = weights[iEdge]*(nodeLast[0]
            - nodeFirst[0]);
       const CFreal dy = weights[iEdge]*(nodeLast[1]
            - nodeFirst[1]);

       CFLogDebugMax( "weights = " << weights[iEdge]
       << "dx = " << dx
       << "dy = " << dy << "\n");

       _l11[firstID] += dx*dx;
       _l12[firstID] += dx*dy;
       _l22[firstID] += dy*dy;
       
       if (!last->isGhost()) {
	 _l11[lastID] += dx*dx;
	 _l12[lastID] += dx*dy;
	 _l22[lastID] += dy*dy;
       }
       
       ++iEdge;
     }
   }
 }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
