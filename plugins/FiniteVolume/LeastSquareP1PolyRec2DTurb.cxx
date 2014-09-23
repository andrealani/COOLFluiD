#include "LeastSquareP1PolyRec2DTurb.hh"
#include "Framework/VolumeIntegrator.hh"
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

MethodStrategyProvider<LeastSquareP1PolyRec2DTurb, CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule> 
leastSquareP1PolyRec2DTurbProvider("LinearLS2DTurb");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DTurb::LeastSquareP1PolyRec2DTurb(const std::string& name) :
  LeastSquareP1PolyRec2D(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DTurb::~LeastSquareP1PolyRec2DTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DTurb::extrapolateImpl(GeometricEntity* const face)
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
  
  for(CFuint iVar = 0; iVar < 6; ++iVar) {
    cf_assert(iVar < _zeroGradient->size());
    
    const CFreal xq = (*coord[0])[XX];
    const CFreal yq = (*coord[0])[YY];
    
    const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ? 
      uX(stateID,iVar,nbEqs)*(xq - stateCoord[0]) + 
      uY(stateID,iVar,nbEqs)*(yq - stateCoord[1]) : 0.0;
    
    const CFreal gradientCoeffNeighbor = (!(*_zeroGradient)[iVar]) ? 
      uX(neighStateIDgrad,iVar,nbEqs)*(xq - neighStateCoord[0]) + 
      uY(neighStateIDgrad,iVar,nbEqs)*(yq - neighStateCoord[1]) : 0.0;
    
    // left reconstructed value (the one inside the current cell)
    (*valuesLR[0])[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
    
    // right reconstructed value (the one inside the neighbor cell)
    (*valuesLR[1])[iVar] = (*neighState)[iVar] + newLimiter[neighStartID + iVar]*gradientCoeffNeighbor;
    
    getBackupValues(0)[iVar] =  (*valuesLR[0])[iVar];
    getBackupValues(1)[iVar] =  (*valuesLR[1])[iVar];
  }

  for(CFuint iVar = 4; iVar < nbEqs; ++iVar) {
      (*valuesLR[0])[iVar] = (*face->getState(LEFT))[iVar];
      (*valuesLR[0])[iVar] = (*face->getState(RIGHT))[iVar];

      getBackupValues(0)[iVar] = (*valuesLR[0])[iVar];
      getBackupValues(1)[iVar] = (*valuesLR[0])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DTurb::extrapolateImpl(GeometricEntity* const face,
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
  
  // number of quadrature points associated with this face
  const vector<Node*>& coord = getCoord();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint startID = stateID*nbEqs;

  const CFreal xq = (*coord[0])[XX];
  const CFreal yq = (*coord[0])[YY];
  copyBackupValues();

  const CFreal gradientCoeffState = uX(stateID,iVar,nbEqs)*
      (xq - stateCoord[0]) + uY(stateID,iVar,nbEqs)*
      (yq - stateCoord[1]);

    // reconstructed value
    if(iVar < 4){
      getValues(leftOrRight)[iVar] = (*state)[iVar] +
        newLimiter[startID + iVar]*gradientCoeffState;
    }
    if( (iVar >=4) && (iVar < nbEqs) ){
      getValues(leftOrRight)[iVar] = (*state)[iVar];
    }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
