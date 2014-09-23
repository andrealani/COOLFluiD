#include "LeastSquareP1PolyRec3DFixYi.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
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

MethodStrategyProvider<LeastSquareP1PolyRec3DFixYi, CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, FiniteVolumeNEQModule> 
leastSquareP1PolyRec3DFixYiProvider("LinearLS3DFixYi");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec3DFixYi::LeastSquareP1PolyRec3DFixYi(const std::string& name) :
  LeastSquareP1PolyRec3D(name),
  m_yiL(),
  m_yiR()
{
}
      
//////////////////////////////////////////////////////////////////////////////
      
LeastSquareP1PolyRec3DFixYi::~LeastSquareP1PolyRec3DFixYi()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3DFixYi::setup()
{
  LeastSquareP1PolyRec3D::setup();
  
  const CFuint nbSpecies = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>()->getNbSpecies();
  
  m_yiL.resize(nbSpecies);
  m_yiR.resize(nbSpecies);
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3DFixYi::extrapolateImpl(GeometricEntity* const face)
{
  FVMCC_PolyRec::baseExtrapolateImpl(face);
  
  // first compute the limiter
  computeFaceLimiter(face);
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  // please note that you work with references !!!
  // don't forget the "&" !!!
  const vector<Node*>& coord = getCoord();
  vector<State*>& valuesLR  = getExtrapolatedValues();
  const State *const state = face->getState(LEFT);
  const CFuint stateID = state->getLocalID();
  const RealVector& stateCoord = state->getCoordinates();
  const State *const neighState = face->getState(RIGHT);
  const RealVector& neighStateCoord = neighState->getCoordinates();
  const CFuint neighStateIDgrad = (!isBoundaryFace()) ? neighState->getLocalID() : stateID;
  
  // number of quadrature points associated with this face
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint startID = stateID*nbEqs;
  const CFuint neighStartID = neighStateIDgrad*nbEqs;
  
  // sanity check
  const CFuint nbSpecies = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>()->getNbSpecies();
  
  CFreal rhoL = 0.0;  
  CFreal rhoR = 0.0;
  for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    const CFreal xq = (*coord[0])[XX];
    const CFreal yq = (*coord[0])[YY];
    const CFreal zq = (*coord[0])[ZZ];
    const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ?
      uX(stateID,iVar,nbEqs)*(xq - stateCoord[XX]) + 
      uY(stateID,iVar,nbEqs)*(yq - stateCoord[YY]) + 
      uZ(stateID,iVar,nbEqs)*(zq - stateCoord[ZZ]) : 0.0;
    
    const CFreal gradientCoeffNeighbor = (!(*_zeroGradient)[iVar]) ?
      uX(neighStateIDgrad,iVar,nbEqs)*(xq - neighStateCoord[XX]) + 
      uY(neighStateIDgrad,iVar,nbEqs)*(yq - neighStateCoord[YY]) + 
      uZ(neighStateIDgrad,iVar,nbEqs)*(zq - neighStateCoord[ZZ]) : 0.0;
    
    // left reconstructed value (the one inside the current cell)
    (*valuesLR[0])[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
    
    // right reconstructed value (the one inside the neighbor cell)
    (*valuesLR[1])[iVar] = (*neighState)[iVar] + newLimiter[neighStartID + iVar]*gradientCoeffNeighbor;
    
    if (iVar < nbSpecies) {
      rhoL += (*valuesLR[0])[iVar];
      rhoR += (*valuesLR[1])[iVar];
      
      m_yiL[iVar] = (*valuesLR[0])[iVar];
      m_yiR[iVar] = (*valuesLR[1])[iVar];
    }
    
    getBackupValues(0)[iVar] = (*valuesLR[0])[iVar];
    getBackupValues(1)[iVar] = (*valuesLR[1])[iVar];
  }
  m_yiL /= rhoL;
  m_yiR /= rhoR;
  
  const CFreal sumYiL = m_yiL.sum();
  const CFreal sumYiR = m_yiR.sum();
  
  // if sum of mass fractions is > 1.0 resort to first order
  if (sumYiL > 1.001) {
    for(CFuint iVar = 0; iVar < nbSpecies; ++iVar) {
      (*valuesLR[0])[iVar] = (*state)[iVar];
      getBackupValues(0)[iVar] = (*valuesLR[0])[iVar];
    }
  }
  
  // if sum of mass fractions is > 1.0 resort to first order
  if (sumYiR > 1.001) {
    for(CFuint iVar = 0; iVar < nbSpecies; ++iVar) {
      (*valuesLR[1])[iVar] = (*neighState)[iVar];
      getBackupValues(1)[iVar] = (*valuesLR[1])[iVar];
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////
      
void LeastSquareP1PolyRec3DFixYi::extrapolateImpl(GeometricEntity* const face,
						  CFuint iVar, CFuint leftOrRight)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
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
  const CFuint startID = stateID*PhysicalModelStack::getActive()->getNbEq();

  const CFreal xq = (*coord[0])[XX];
  const CFreal yq = (*coord[0])[YY];
  const CFreal zq = (*coord[0])[ZZ];
  copyBackupValues();
  
  const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ?
    uX(stateID,iVar,nbEqs)*(xq - stateCoord[XX]) + 
    uY(stateID,iVar,nbEqs)*(yq - stateCoord[YY]) + 
    uZ(stateID,iVar,nbEqs)*(zq - stateCoord[ZZ]) : 0.0;
  
  // reconstructed value
  getValues(leftOrRight)[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

 } // namespace COOLFluiD
