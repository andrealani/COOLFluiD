#include "LeastSquareP1PolyRec2DBcFix.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Common/NotImplementedException.hh"

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

MethodStrategyProvider<LeastSquareP1PolyRec2DBcFix, 
		       CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeModule> 
leastSquareP1PolyRec2DBcFixProvider("LinearLS2DBcFix");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DBcFix::LeastSquareP1PolyRec2DBcFix
(const std::string& name) :
  LeastSquareP1PolyRec2D(name),
  m_drXiXg(0.),
  m_drXqXg(0.),
  m_drXqXi(0.)
{
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DBcFix::~LeastSquareP1PolyRec2DBcFix()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DBcFix::extrapolateImpl(GeometricEntity* const face)
{
  if (!isBoundaryFace()){
    LeastSquareP1PolyRec2D::extrapolateImpl(face);
  }
  else {
    cf_assert(isBoundaryFace());
    DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
    
    // please note that you work with references !!!
    // don't forget the "&" !!!
    const vector<Node*>& coord = getCoord();
    vector<State*>& valuesLR  = getExtrapolatedValues();
    const State *const state = face->getState(LEFT);
    const State *const neighState = face->getState(RIGHT);
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    cf_assert(_zeroGradient != CFNULL);
    
    // m_dr* could be put in the FluxSplitterData and computed by 
    // the BC condition in some cases
    // distance between internal and ghost state 
    m_drXiXg = MathFunctions::getDistance(state->getCoordinates(), neighState->getCoordinates());
    
    // we suppose to have just one quadrature point
    m_drXqXg = MathFunctions::getDistance(*coord[0], neighState->getCoordinates());
    
    // we suppose to have just one quadrature point
    m_drXqXi = MathFunctions::getDistance(*coord[0], state->getCoordinates());

    for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ? 
	((*state)[iVar] - (*neighState)[iVar])/m_drXiXg : 0.0;
      
      // no limiter for the moment ...

      // left reconstructed value (the one inside the current cell)
      (*valuesLR[0])[iVar] = (*state)[iVar] + gradientCoeffState*m_drXqXi;

     //  newLimiter[startID + iVar]*gradientCoeffState*m_drXqXi;
     // (*leftValues[0])[iVar] = (*state)[iVar] + gradientCoeffState;
      
      // right reconstructed value (the one inside the neighbor cell)
      // (*rightValues[0])[iVar] = (*neighState)[iVar] - gradientCoeffState;
      (*valuesLR[1])[iVar] = (*neighState)[iVar] + gradientCoeffState*m_drXqXg;
      
      //  newLimiter[startID + iVar]*gradientCoeffState*m_drXqXg;
      
      getBackupValues(0)[iVar] =  (*valuesLR[0])[iVar];
      getBackupValues(1)[iVar] =  (*valuesLR[1])[iVar];
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DBcFix::extrapolateImpl
(GeometricEntity* const face, CFuint iVar, CFuint leftOrRight)
{
  /// RECHECK THIS !!!!!
  if (!isBoundaryFace()){
    LeastSquareP1PolyRec2D::extrapolateImpl(face, iVar, leftOrRight);
  }
  else {
    DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
    
    const State *const state = face->getState(leftOrRight);
    CFuint stateID = face->getState(LEFT)->getLocalID();
    
    const CFuint otherLR = (leftOrRight == 0) ? 1 : 0;
    const State *const neighState = face->getState(otherLR);
    const CFuint startID = stateID*PhysicalModelStack::getActive()->getNbEq();
    const vector<Node*>& coord = getCoord();

    // m_dr* could be put in the FluxSplitterData and computed by 
    // the BC condition in some cases
    // distance between internal and ghost state 
    m_drXiXg = MathFunctions::getDistance(state->getCoordinates(),neighState->getCoordinates());
    
    // we suppose to have just one quadrature point
    m_drXqXg = MathFunctions::getDistance(*coord[0], neighState->getCoordinates());
    
    copyBackupValues();
    
    const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ? 
      ((*state)[iVar] - (*neighState)[iVar])/m_drXiXg*m_drXqXg : 0.0;
    
    // reconstructed value
    getValues(leftOrRight)[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
