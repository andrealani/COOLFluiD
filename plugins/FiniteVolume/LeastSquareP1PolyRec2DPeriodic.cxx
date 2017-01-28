#include "FiniteVolume/LeastSquareP1PolyRec2DPeriodic.hh"
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

MethodStrategyProvider<LeastSquareP1PolyRec2DPeriodic, CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeModule> 
leastSquareP1PolyRec2DPeriodicProvider("LinearLS2DPeriodic");

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >
    ("PeriodicTRS","Names of the periodic TRS(s).");
}
      
//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DPeriodic::LeastSquareP1PolyRec2DPeriodic(const std::string& name) :
  LeastSquareP1PolyRec2D(name),
  m_isPeriodicFace()
{
  addConfigOptionsTo(this);
  
  m_periodicTRS = std::vector<std::string>();
  setParameter("PeriodicTRS", &m_periodicTRS);
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DPeriodic::~LeastSquareP1PolyRec2DPeriodic()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LeastSquareP1PolyRec2DPeriodic::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    LeastSquareP1PolyRec2D::needsSockets();
  
  // add something?
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::computeGradients()
{
  LeastSquareP1PolyRec2D::computeGradients();
  
  // communicate the gradients corresponding to the cells on the 
  // matching side of local periodic faces and store them in a local map or array
  // this has to be done in a smart way with minimal memory allocation and 
  // MPI-friendly 
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::extrapolateImpl(GeometricEntity* const face)
{
  if (!m_isPeriodicFace[face->getID()]) { 
    LeastSquareP1PolyRec2D::extrapolateImpl(face);
  }
  else {
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
      (*valuesLR[0])[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
      
      // right reconstructed value (the one inside the neighbor cell)
      (*valuesLR[1])[iVar] = (*neighState)[iVar] + newLimiter[neighStartID + iVar]*gradientCoeffNeighbor;
      
      getBackupValues(0)[iVar] =  (*valuesLR[0])[iVar];
      getBackupValues(1)[iVar] =  (*valuesLR[1])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::extrapolateImpl(GeometricEntity* const face,
					     CFuint iVar, CFuint leftOrRight)
{
  if (!m_isPeriodicFace[face->getID()]) {
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
    getValues(leftOrRight)[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
  }
  else {
    LeastSquareP1PolyRec2DPeriodic::extrapolateImpl(face, iVar, leftOrRight);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::setup()
{
  LeastSquareP1PolyRec2D::setup();
  
  // AL: consider adding attributes to TopologicalRegionSet 
  // (e.g. isPeriodic, isWall, isPartition) in order to avoid to always
  // create flag arrays like m_isPeriodicFace
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  m_isPeriodicFace.resize(nbFaces, false);
  
  /// here flag the periodic faces
  for (CFuint iTRS = 0; iTRS < m_periodicTRS.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = 
      MeshDataStack::getActive()->getTrs(m_periodicTRS[iTRS]);
    const CFuint nbLocalFaces = trs->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbLocalFaces; ++iFace) {
      m_isPeriodicFace[trs->getLocalGeoID(iFace)] = true;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
