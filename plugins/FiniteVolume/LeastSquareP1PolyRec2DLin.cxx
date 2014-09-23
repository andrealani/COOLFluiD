#include "LeastSquareP1PolyRec2DLin.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
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

MethodStrategyProvider<LeastSquareP1PolyRec2DLin, CellCenterFVMData, PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule> LeastSquareP1PolyRec2DLinProvider("LinearLS2DLin");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DLin::LeastSquareP1PolyRec2DLin(const std::string& name) :
  FVMCC_PolyRecLin(name),
  socket_stencil("stencil"),
  socket_weights("weights"),
  socket_uX("uX"),
  socket_uY("uY"),
  socket_linearizedStates("linearizedStates"),
  socket_linearizedGhostStates("linearizedGhostStates"),
  _l11(),
  _l12(),
  _l22(),
  _lf1(),
  _lf2()
{
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DLin::~LeastSquareP1PolyRec2DLin()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DLin::configure ( Config::ConfigArgs& args )
{
  FVMCC_PolyRecLin::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LeastSquareP1PolyRec2DLin::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    FVMCC_PolyRecLin::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_stencil);
  result.push_back(&socket_weights);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_linearizedStates);
  result.push_back(&socket_linearizedGhostStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DLin::computeGradients()
{
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
      const State* const first = states[iState];
      const CFuint stencilSize = stencil[iState].size();
      // loop over the neighbor cells belonging to the chosen stencil
      for(CFuint in = 0; in < stencilSize; ++in) {
	const State* const last = stencil[iState][in];
	const CFuint lastID = last->getLocalID();
	const CFuint firstID = first->getLocalID();
	cf_assert(firstID != lastID);

	if (lastID > firstID) {
	  // consider the next edge
	  const RealVector& nodeFirst = first->getCoordinates();
	  const RealVector& nodeLast = last->getCoordinates();

	  const CFreal dx = weights[iEdge]*(nodeLast[0]
					     - nodeFirst[0]);
	  const CFreal dy = weights[iEdge]*(nodeLast[1]
					     - nodeFirst[1]);

	  const CFreal du = weights[iEdge]*((*last)[iVar] - (*first)[iVar]);
	  CFLogDebugMax( "du = " << du << "\n");

	  _lf1[firstID] += dx*du;
	  _lf1[lastID]  += dx*du;
	  _lf2[firstID] += dy*du;
	  _lf2[lastID]  += dy*du;

	  ++iEdge;
	}
      }
    }

    for(CFuint iState = 0; iState < nbStates; ++iState) {
      const CFreal det = _l11[iState]*_l22[iState] - _l12[iState]*_l12[iState];
      uX(iState,iVar,nbEquations) = (_l22[iState]*_lf1[iState] -
				     _l12[iState]*_lf2[iState])/det;
      uY(iState,iVar,nbEquations) = (_l11[iState]*_lf2[iState] -
			  _l12[iState]*_lf1[iState])/det;

      CFLogDebugMax( "det = " << det
		     << ", l11 = " << _l11[iState]
		     << ", l12 = " << _l12[iState]
		     << ", l22 = " << _l22[iState]
		     << ", lf1 = " << _lf1[iState]
		     << ", lf2 = " << _lf2[iState]
		     << ", uX = " << uX(iState,iVar,nbEquations)
		     << ", uY = " << uY(iState,iVar,nbEquations)
		     << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DLin::extrapolateImpl(GeometricEntity* const face)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  // please note that you work with references !!!
  // don't forget the "&" !!!
  const vector<Node*>& coord = getCoord();
  vector<State*>& leftValues  = getValues(LEFT);
  vector<State*>& rightValues = getValues(RIGHT);
  vector<RealVector>& backUpLeftValues = getBackupValues(LEFT);
  vector<RealVector>& backUpRightValues = getBackupValues(RIGHT);

  VolumeIntegratorImpl* const impl =
    getMethodData().getVolumeIntegrator()->getSolutionIntegrator(face);

  const std::valarray<CFreal>& coeff = impl->getCoeff();

  const State *const state = face->getState(LEFT);
  const CFuint stateID = state->getLocalID();
  const RealVector& stateCoord = state->getCoordinates();

  const State *const neighState = face->getState(RIGHT);
  const CFuint neighStateIDgrad = (!isBoundaryFace()) ?
    neighState->getLocalID() : stateID;
  const RealVector& neighStateCoord = neighState->getCoordinates();

  // number of quadrature points associated with this face
  const CFuint nbPoints = coeff.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint startID = stateID*nbEqs;
  const CFuint neighStartID = neighStateIDgrad*nbEqs;
  
  for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFuint ir = 0;

    for (CFuint ip = 0; ip < nbPoints; ++ip, ++ir) {
      const CFreal xq = (*coord[ir])[0];
      const CFreal yq = (*coord[ir])[1];
      const CFreal gradientCoeffState = uX(stateID,iVar,nbEqs)*
	(xq - stateCoord[0]) + uY(stateID,iVar,nbEqs)*
	(yq - stateCoord[1]);

      const CFreal gradientCoeffNeighbor = uX(neighStateIDgrad,iVar,nbEqs)*
	(xq - neighStateCoord[0]) + uY(neighStateIDgrad,iVar,nbEqs)*
	(yq - neighStateCoord[1]);
      
      // left reconstructed value (the one inside the current cell)
      (*leftValues[ip])[iVar] = (*state)[iVar] +
	newLimiter[startID + iVar]*gradientCoeffState;
      
      // right reconstructed value (the one inside the neighbor cell)
      (*rightValues[ip])[iVar] = (*neighState)[iVar] +
	newLimiter[neighStartID + iVar]*gradientCoeffNeighbor;
      
      backUpLeftValues[ip][iVar] = (*leftValues[ip])[iVar];
      backUpRightValues[ip][iVar] = (*rightValues[ip])[iVar];
    }
  }


  // for each quadrature point, set
  // left state  = current state
  // right state = neighbor state
  vector<State*>& leftValues2 = getValues2(LEFT);
  vector<State*>& rightValues2 = getValues2(RIGHT);
  vector<RealVector>& backUpLeftValues2 = getBackupValues2(LEFT);
  vector<RealVector>& backUpRightValues2 = getBackupValues2(RIGHT);

  const bool isGhostLeft  = (*face->getState(LEFT)).isGhost();
  const bool isGhostRight = (*face->getState(RIGHT)).isGhost();

  State* leftState2(CFNULL);
  State* rightState2(CFNULL);

  if(isGhostLeft) leftState2 = (socket_linearizedGhostStates.getDataHandle())[(*face->getState(LEFT)).getLocalID()];
  else leftState2 = (socket_linearizedStates.getDataHandle())[(*face->getState(LEFT)).getLocalID()];

  if(isGhostRight) rightState2 = (socket_linearizedGhostStates.getDataHandle())[(*face->getState(RIGHT)).getLocalID()];
  else rightState2 = (socket_linearizedStates.getDataHandle())[(*face->getState(RIGHT)).getLocalID()];


  for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFuint ir = 0;

    for (CFuint ip = 0; ip < nbPoints; ++ip, ++ir) {
      const CFreal xq = (*coord[ir])[0];
      const CFreal yq = (*coord[ir])[1];
      const CFreal gradientCoeffState = uX(stateID,iVar,nbEqs)*
	(xq - stateCoord[0]) + uY(stateID,iVar,nbEqs)*
	(yq - stateCoord[1]);

      const CFreal gradientCoeffNeighbor = uX(neighStateIDgrad,iVar,nbEqs)*
	(xq - neighStateCoord[0]) + uY(neighStateIDgrad,iVar,nbEqs)*
	(yq - neighStateCoord[1]);

      // left reconstructed value (the one inside the current cell)
      (*leftValues2[ip])[iVar] = (*leftState2)[iVar] +
	newLimiter[startID + iVar]*gradientCoeffState;

      // right reconstructed value (the one inside the neighbor cell)
      (*rightValues2[ip])[iVar] = (*rightState2)[iVar] +
        newLimiter[neighStartID + iVar]*gradientCoeffNeighbor;

      backUpLeftValues2[ip][iVar] = (*leftValues2[ip])[iVar];
      backUpRightValues2[ip][iVar] = (*rightValues2[ip])[iVar];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DLin::extrapolateImpl(GeometricEntity* const face,
					     CFuint iVar, CFuint leftOrRight)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  const State *const state = face->getState(leftOrRight);
  CFuint stateID = 0;
  if (leftOrRight == LEFT) {
    stateID = state->getLocalID();
  }
  else {
    stateID = (!isBoundaryFace()) ? state->getLocalID() :
      face->getState(LEFT)->getLocalID();
  }
  const RealVector& stateCoord = state->getCoordinates();

  VolumeIntegratorImpl* const impl =
   getMethodData().getVolumeIntegrator()->getSolutionIntegrator(face);

  const std::valarray<CFreal>& coeff = impl->getCoeff();
  // number of quadrature points associated with this face
  const CFuint nbPoints = coeff.size();
  const vector<Node*>& coord = getCoord();
  vector<State*>& rvalues  = getValues(leftOrRight);
  const CFuint startID = stateID*PhysicalModelStack::getActive()->getNbEq();
    
  CFuint ir = 0;
  for (CFuint ip = 0; ip < nbPoints; ++ip, ++ir) {
    copyBackupValues(ip);

    const CFreal xq = (*coord[ir])[0];
    const CFreal yq = (*coord[ir])[1];
    const CFreal gradientCoeffState = uX(stateID,iVar,nbEqs)*
      (xq - stateCoord[0]) + uY(stateID,iVar,nbEqs)*
      (yq - stateCoord[1]);

    // reconstructed value
    (*rvalues[ip])[iVar] = (*state)[iVar] +
      newLimiter[startID + iVar]*gradientCoeffState;
  }



  const bool isGhost = state->isGhost();

  State* state2(CFNULL);

  if(isGhost) state2 = (socket_linearizedGhostStates.getDataHandle())[(*face->getState(leftOrRight)).getLocalID()];
  else state2 = (socket_linearizedStates.getDataHandle())[(*face->getState(leftOrRight)).getLocalID()];

  // for each quadrature point, set
  // reconstructed variable = same variable in the current state
  vector<State*>& values2 = getValues2(leftOrRight);

  ir = 0;
  for (CFuint ip = 0; ip < nbPoints; ++ip, ++ir) {
    copyBackupValues2(ip);

    const CFreal xq = (*coord[ir])[0];
    const CFreal yq = (*coord[ir])[1];
    const CFreal gradientCoeffState = uX(stateID,iVar,nbEqs)*
      (xq - stateCoord[0]) + uY(stateID,iVar,nbEqs)*
      (yq - stateCoord[1]);
    
    // reconstructed value
    (*values2[ip])[iVar] = (*state2)[iVar] +
      newLimiter[startID + iVar]*gradientCoeffState;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DLin::updateWeights()
{

  FVMCC_PolyRecLin::updateWeights();

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
     const CFuint lastID = last->getLocalID();
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
       _l11[lastID] += dx*dx;
       _l12[lastID] += dx*dy;
       _l22[lastID] += dy*dy;

       ++iEdge;
     }
   }
 }

}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DLin::setup()
{
  FVMCC_PolyRecLin::setup();

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
     const CFuint lastID = last->getLocalID();
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
       _l11[lastID] += dx*dx;
       _l12[lastID] += dx*dy;
       _l22[lastID] += dy*dy;

       ++iEdge;
     }
   }
 }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
