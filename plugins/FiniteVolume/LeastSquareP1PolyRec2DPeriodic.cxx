#include "FiniteVolume/LeastSquareP1PolyRec2DPeriodic.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_BCPeriodic.hh"

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
    ("PeriodicBCNames","Names of the periodic BC commands.");
}
      
//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec2DPeriodic::LeastSquareP1PolyRec2DPeriodic(const std::string& name) :
  LeastSquareP1PolyRec2D(name),
  m_mapGeoToTrs(CFNULL),
  m_isPeriodicFace(),
  m_periodicGradients(),
  m_mapTrs2BC()
{
  addConfigOptionsTo(this);
  
  m_periodicBCNames = std::vector<std::string>();
  setParameter("PeriodicBCNames", &m_periodicBCNames);
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

  CFuint counter = 0;    
  SafePtr<CFMap<CFuint, FVMCC_BC*> > bcMap = getMethodData().getMapBC();
  for (CFuint iTRS = 0; iTRS < bcMap->size(); ++iTRS) {
    if (binary_search(m_periodicBCNames.begin(), 
		      m_periodicBCNames.end(),
		      (*bcMap)[iTRS]->getName())) {
      (*bcMap)[iTRS]->transferGradientsData();
      counter++;
    }
  }
  cf_assert(m_periodicBCNames.size() == counter);
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
    const RealVector& neighStateCoord = neighState->getCoordinates();

    // number of quadrature points associated with this face
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint startID = stateID*nbEqs;
    cf_assert(_zeroGradient != CFNULL);
    
    const CFreal xq = (*coord[0])[XX];
    const CFreal yq = (*coord[0])[YY];
    
    if (_vars.size() > 1) {
      _vFunction.evaluate((*coord[0]),_gradientCoeff);
    }
    
    CFreal* periodicLimiter = CFNULL;
    TopologicalRegionSet *const trs = &*m_mapGeoToTrs->getTrs(face->getID());
    m_mapTrs2BC.find(trs)->second->computePeriodicGradient(face, m_periodicGradients, periodicLimiter);
    
    for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      cf_assert(iVar < _zeroGradient->size());
      
      // AL: gradientCoeff prevails on the interactive gradientFactor
      if (_vars.size() > 1) {_gradientFactor = _gradientCoeff[iVar];}
      
      const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ? 
	_gradientFactor*(uX(stateID,iVar,nbEqs)*(xq - stateCoord[0]) + 
			 uY(stateID,iVar,nbEqs)*(yq - stateCoord[1])) : 0.0;
      
      const CFreal gradientCoeffNeighbor = (!(*_zeroGradient)[iVar]) ? 
	_gradientFactor*(m_periodicGradients[XX][iVar]*(xq - neighStateCoord[0]) + 
			 m_periodicGradients[YY][iVar]*(yq - neighStateCoord[1])) : 0.0;
      
      // left reconstructed value (the one inside the current cell)
      (*valuesLR[0])[iVar] = (*state)[iVar] + newLimiter[startID + iVar]*gradientCoeffState;
      
      // right reconstructed value (the one inside the neighbor cell)
      (*valuesLR[1])[iVar] = (*neighState)[iVar] + periodicLimiter[iVar]*gradientCoeffNeighbor;
      
      getBackupValues(0)[iVar] =  (*valuesLR[0])[iVar];
      getBackupValues(1)[iVar] =  (*valuesLR[1])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::extrapolateImpl(GeometricEntity* const face,
						     CFuint iVar, CFuint leftOrRight)
{
  CFLog(DEBUG_MED, "LeastSquareP1PolyRec2DPeriodic::extrapolateImpl() => start\n");
  
  if (!m_isPeriodicFace[face->getID()] ||
      (m_isPeriodicFace[face->getID()] && leftOrRight == LEFT))  {
    LeastSquareP1PolyRec2D::extrapolateImpl(face, iVar, leftOrRight);
    CFLog(DEBUG_MED, "LeastSquareP1PolyRec2DPeriodic::extrapolateImpl() => end\n");
  }
  else {
    cf_assert(m_isPeriodicFace[face->getID()] && leftOrRight == RIGHT);
    DataHandle<CFreal> uX = socket_uX.getDataHandle();
    DataHandle<CFreal> uY = socket_uY.getDataHandle();
    DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
    const State *const state = face->getState(leftOrRight);
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint stateID = face->getState(LEFT)->getLocalID();
    const RealVector& stateCoord = state->getCoordinates();
    
    // number of quadrature points associated with this face
    const vector<Node*>& coord = getCoord();
    const CFreal xq = (*coord[0])[XX];
    const CFreal yq = (*coord[0])[YY];
    copyBackupValues();
    
    if (_vars.size() > 1) {
      _vFunction.evaluate(iVar,(*coord[0]), _gradientFactor);
    }
    
    CFreal* periodicLimiter = CFNULL;
    TopologicalRegionSet *const trs = &*m_mapGeoToTrs->getTrs(face->getID());
    m_mapTrs2BC.find(trs)->second->computePeriodicGradient(face, m_periodicGradients, periodicLimiter);
    
    // AL: gradientCoeff prevails on the interactive gradientFactor
    const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ?
      _gradientFactor*(m_periodicGradients[XX][iVar]*(xq - stateCoord[XX]) + 
		       m_periodicGradients[YY][iVar]*(yq - stateCoord[YY])) : 0.0;
    
    // reconstructed value
    getValues(leftOrRight)[iVar] = (*state)[iVar] + periodicLimiter[iVar]*gradientCoeffState;
  }
  
  CFLog(DEBUG_MED, "LeastSquareP1PolyRec2DPeriodic::extrapolateImpl() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec2DPeriodic::setup()
{
  LeastSquareP1PolyRec2D::setup();
  
  m_mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  cf_assert(m_mapGeoToTrs != CFNULL);
  
  sort(m_periodicBCNames.begin(), m_periodicBCNames.end());
  
  m_periodicGradients.resize(PhysicalModelStack::getActive()->getDim());
  
  // AL: consider adding attributes to TopologicalRegionSet 
  // (e.g. isPeriodic, isWall, isPartition) in order to avoid to always
  // create flag arrays like m_isPeriodicFace
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  m_isPeriodicFace.resize(nbFaces, false);
  
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr<CFMap<CFuint, FVMCC_BC*> > bcMap = getMethodData().getMapBC();
  
  CFuint counter = 0;
  for (CFuint iTRS = 0; iTRS < bcMap->size(); ++iTRS) {
    if (binary_search(m_periodicBCNames.begin(), 
		      m_periodicBCNames.end(),
		      (*bcMap)[iTRS]->getName())) {
      const CFuint trsID = bcMap->getKey(iTRS);
      m_mapTrs2BC[&*trs[trsID]] = dynamic_cast<BCPeriodic*>((*bcMap)[iTRS]);
      const CFuint nbLocalFaces = trs[trsID]->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbLocalFaces; ++iFace) {
	const CFuint faceID = trs[trsID]->getLocalGeoID(iFace);
	m_isPeriodicFace[faceID] = true;
      }
      counter++;
    }
  }
  
  cf_assert(m_periodicBCNames.size() == counter);
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
