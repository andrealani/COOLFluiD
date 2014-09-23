#include "LaxFriedFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LaxFriedFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
stdLaxFriedFluxProvider("LaxFried");

//////////////////////////////////////////////////////////////////////////////

void LaxFriedFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >
    ("DiffCoeff", "Diffusion reduction coefficient");
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedFlux::LaxFriedFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _solutionStates(CFNULL),
  _statesLR(2),
  _sumFlux(),
  _rightEv(),
  _leftEv(),
  _tempUnitNormal()
{
  addConfigOptionsTo(this);
  _currentDiffRedCoeff = 1.0;
  setParameter("DiffCoeff", &_currentDiffRedCoeff); 
}
      
//////////////////////////////////////////////////////////////////////////////
      
void LaxFriedFlux::configure ( Config::ConfigArgs& args )
{
  FVMCC_FluxSplitter::configure(args);
}

/////////////////////////////////////////////////////////////////////////////

LaxFriedFlux::~LaxFriedFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedFlux::setup()
{
  FVMCC_FluxSplitter::setup();
  
  _sumFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedFlux::compute(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  CellCenterFVMData& data = this->getMethodData(); 
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  
  // flux for the right and left state
  _sumFlux = updateVarSet->getFlux()(pdata[1], unitNormal);
  _sumFlux += updateVarSet->getFlux()(pdata[0], unitNormal);
    
  updateVarSet->computeEigenValues(pdata[1], unitNormal, _rightEv);
  updateVarSet->computeEigenValues(pdata[0], unitNormal, _leftEv);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFreal aR = 0.0;
  CFreal aL = 0.0;
  for (CFuint i = 0; i < nbEqs; ++i) {
    aR = max(aR, std::abs(_rightEv[i]));
    aL = max(aL, std::abs(_leftEv[i]));
  }
  const CFreal a = max(aR,aL);
  
  _statesLR[0] = &polyRec->getCurrLeftState();
  _statesLR[1] = &polyRec->getCurrRightState();
  cf_assert(*_statesLR[0] == polyRec->getCurrLeftState());
  cf_assert(*_statesLR[1] == polyRec->getCurrRightState());
  
  if (!getMethodData().reconstructSolVars()) {
    _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  }
  else {
    _solutionStates = &_statesLR;
  }
  
  // you must work with references (no copying allowed) !!!!
  State& leftState  = *(*_solutionStates)[0];
  State& rightState = *(*_solutionStates)[1];
  const CFreal aDiff = a*getReductionCoeff();
 
  result = 0.5*(_sumFlux - aDiff*(rightState - leftState));

  // compute update coefficient
  if (!getMethodData().isPerturb()) {
    CellCenterFVMData& data = this->getMethodData(); 
    GeometricEntity& face = *data.getCurrentFace();
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      polyRec->nbQPoints();
    
    // left contribution to update coefficient
    // CFreal maxEV = updateVarSet->getMaxEigenValue(pdata[0], unitNormal);
    
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(_leftEv.max(), 0.0)*faceArea;
    // updateCoeff[leftID] += max(maxEV, 0.)*faceArea;
    
    /*if (leftID == 0) {    
      const CFuint stateR = getMethodData().getCurrentFace()->getState(1)->getLocalID();
      std::cout << endl << "stateIDs [" << leftID << "," << stateR << "]" << std::endl;
      // std::cout.precision(12);std::cout << "centerNodes[0]    = " << getMethodData().getCurrentFace()->getState(0)->getCoordinates() << endl;
      std::cout.precision(12);std::cout << "dataR   = " << pdata[1] << std::endl;
      // std::cout.precision(12);std::cout << "eigenR  = " << _rightEv << std::endl;
      std::cout.precision(12);std::cout << "normalR = " << unitNormal << std::endl;
      std::cout.precision(12);std::cout << "dataL   = " << pdata[0] << std::endl;
      // std::cout.precision(12);std::cout << "eigenL  = " << _leftEv << std::endl;
      std::cout.precision(12);std::cout << "normalL = " << unitNormal << std::endl;
      std::cout.precision(12);std::cout << "lambda = " << max(_leftEv.max(), 0.0) << ", area =" << faceArea << std::endl;
      std::cout.precision(12);std::cout << "updateCoeff [" <<leftID << "] = " << max(_leftEv.max(), 0.0)*faceArea << std::endl;
      std::cout.precision(12);;std::cout << "UR       = " << rightState << std::endl;
      std::cout.precision(12);;std::cout << "UL       = " << leftState << std::endl;
      std::cout.precision(12);;std::cout << "_sumFlux = " << _sumFlux << std::endl;
      std::cout.precision(12);;std::cout << "aDiff    = " << aDiff << std::endl;
     }*/
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      
      _tempUnitNormal = -1.0*unitNormal;
      const CFreal maxEV = updateVarSet->getMaxEigenValue(pdata[1],_tempUnitNormal);
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, 0.)*faceArea;
      
      ////
      /*if (rightID == 0) {    
	const CFuint stateR = getMethodData().getCurrentFace()->getState(0)->getLocalID();
	std::cout << endl << "stateIDs [" << rightID << "," << stateR << "]" << std::endl;
	// std::cout.precision(12);std::cout << "centerNodes[0]    = " << getMethodData().getCurrentFace()->getState(1)->getCoordinates() << endl;
	std::cout.precision(12);std::cout << "dataR   = " << pdata[0] << std::endl;
	// std::cout.precision(12);std::cout << "eigenR  = " << _leftEv << std::endl;
	std::cout.precision(12);std::cout << "normalR = " << unitNormal << std::endl;
	std::cout.precision(12);std::cout << "dataL   = " << pdata[1] << std::endl;
	// std::cout.precision(12);std::cout << "eigenL  = " << _rightEv << std::endl;
	std::cout.precision(12);std::cout << "normalL = " << _tempUnitNormal << std::endl;
	std::cout.precision(12);std::cout << "lambda = " << maxEV << ", area =" << faceArea << std::endl;
	std::cout.precision(12);std::cout << "updateCoeff [" <<rightID << "] = " << max(maxEV, 0.)*faceArea << std::endl;
	std::cout.precision(12);;std::cout << "UR       = " << leftState << std::endl;
	std::cout.precision(12);;std::cout << "UL       = " << rightState << std::endl;
	std::cout.precision(12);;std::cout << "_sumFlux = " << _sumFlux << std::endl;
	std::cout.precision(12);;std::cout << "aDiff    = " << aDiff << std::endl;
	}*/
      
      
      
    }
  } 
}

//////////////////////////////////////////////////////////////////////////////

CFreal LaxFriedFlux::getReductionCoeff()
{
  _currentDiffRedCoeff = std::min(_currentDiffRedCoeff, getDissipationControlCoeff());
  return _currentDiffRedCoeff;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
