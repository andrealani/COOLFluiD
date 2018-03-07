#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

void FVMCC_FluxSplitter::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< std::string >("DiffCoeffDef","Definition of the dissipation control function."); 
  options.addConfigOption< std::vector<std::string> >("DiffCoeffVars","Input variables for the dissipation control function.");
}
      
//////////////////////////////////////////////////////////////////////////////

FVMCC_FluxSplitter::FVMCC_FluxSplitter(const std::string& name) :
  FluxSplitter<CellCenterFVMData>(name),
  _dissipationControlCoeff(1.),
  socket_normals("normals"),  
  socket_faceAreas("faceAreas"),
  socket_updateCoeff("updateCoeff"),
  socket_isOutward("isOutward"),
  socket_nstates("nstates"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nodes("nodes"), 
  m_fluxDerivative(),
  _faceArea(0.),
  _firstResidual(0.0),
  _lastResidual(0.0),
  _maxResidual(0.0),
  _tmpJacobMatrix(),
  _dissipationControlInput(),
  _dissipationControlFunction()
{
  addConfigOptionsTo(this);
  
  _dissipationControlDef = "";
  setParameter("DiffCoeffDef",&_dissipationControlDef);
  
  _dissipationControlVars = vector<string>();
  setParameter("DiffCoeffVars",&_dissipationControlVars);
}
  
//////////////////////////////////////////////////////////////////////////////

FVMCC_FluxSplitter::~FVMCC_FluxSplitter()
{
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_FluxSplitter::setup()
{
  FluxSplitter<CellCenterFVMData>::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _tmpJacobMatrix.resize(nbEqs, nbEqs, 0.);
  
  if (_dissipationControlDef != "") {
    cf_assert(_dissipationControlVars.size() > 0);
    // i,r,ri,rl,rmax,cfl
    _dissipationControlInput.resize(_dissipationControlVars.size());
  }
  
  /// build mapping to member functions computing partial derivatives
  // this list should include all variables with respect to which partial derivatives are needed
  // this could be made dynamic using functor objects
  map<string, FluxDerivative> mapVars;
  mapVars["p"]      = &FVMCC_FluxSplitter::dFdP;
  mapVars["u"]      = &FVMCC_FluxSplitter::dFdU;
  mapVars["v"]      = &FVMCC_FluxSplitter::dFdV;
  mapVars["w"]      = &FVMCC_FluxSplitter::dFdW;
  mapVars["T"]      = &FVMCC_FluxSplitter::dFdT;
  mapVars["rho"]    = &FVMCC_FluxSplitter::dFdRho;
  mapVars["K"]      = &FVMCC_FluxSplitter::dFdK;  
  mapVars["Omega"]  = &FVMCC_FluxSplitter::dFdK; 
  mapVars["Nuitil"] = &FVMCC_FluxSplitter::dFdK;
  mapVars["rhoU"]   = &FVMCC_FluxSplitter::dFdRhoU;
  mapVars["rhoV"]   = &FVMCC_FluxSplitter::dFdRhoV;
  mapVars["rhoW"]   = &FVMCC_FluxSplitter::dFdRhoW;
  mapVars["rhoE"]   = &FVMCC_FluxSplitter::dFdRhoE;
  
  // here functions computing flux derivatives are stored in the right order, ready to be applied
  const vector<string>& updateVars = this->getMethodData().getUpdateVar()->getVarNames();
  for (CFuint i = 0; i < updateVars.size(); ++i) {
    if (mapVars.count(updateVars[i]) > 0) {
      m_fluxDerivative.push_back(mapVars.find(updateVars[i])->second);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_FluxSplitter::computeFlux(RealVector& result)
{
  SafePtr<GeometricEntity> currFace = getMethodData().getCurrentFace();
  const CFuint faceID = currFace->getID();
  _faceArea = socket_faceAreas.getDataHandle()[faceID];
  
  // integrate the flux computed at the quadrature points
  if (!getMethodData().useAnalyticalConvJacob()) {
    integrateFluxOnly(result);
  } 
  else {
    const bool isBFace =  currFace->getState(1)->isGhost();
      (!isBFace) ? integrateFluxAndJacob(result) : integrateFluxOnly(result);
  }
}
      
//////////////////////////////////////////////////////////////////////////////
 
void FVMCC_FluxSplitter::configure ( Config::ConfigArgs& args )
{
  FluxSplitter<CellCenterFVMData>::configure(args);
  
  if (_dissipationControlDef != "") {
    if (_dissipationControlVars.size() == 0){
      _dissipationControlVars.resize(6);
      _dissipationControlVars[0] = "i";
      _dissipationControlVars[1] = "r";
      _dissipationControlVars[2] = "ri";
      _dissipationControlVars[3] = "rl";
      _dissipationControlVars[4] = "rmax";
      _dissipationControlVars[5] = "cfl";
    }
    
    vector<string> functions(1, _dissipationControlDef);
    _dissipationControlFunction.setFunctions(functions);
    _dissipationControlFunction.setVariables(_dissipationControlVars);
    try {
      _dissipationControlFunction.parse();
    }
    catch (Common::ParserException& e) {
      CFLog(WARN, e.what() << "\n");
      throw; // rethrow the exception to signal the error to the user
    }
  }
}
 
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FVMCC_FluxSplitter::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result  = FluxSplitter<CellCenterFVMData>::needsSockets();
  
  result.push_back(&socket_normals);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_nstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_FluxSplitter::evaluateDissipationControlFunction()
{  
  // store the first residual
  const CFreal currIter = SubSystemStatusStack::getActive()->getNbIter();
  const CFreal currResidual = SubSystemStatusStack::getActive()->getResidual();
  if (currIter == 1) {
    _firstResidual = SubSystemStatusStack::getActive()->getResidual();
  }
  
  cf_assert(_dissipationControlInput.size() == 6);
  _dissipationControlInput[0] = currIter;
  _dissipationControlInput[1] = currResidual;
  _dissipationControlInput[2] = _firstResidual;
  _dissipationControlInput[3] = _lastResidual;
  _dissipationControlInput[4] = _maxResidual;
  _dissipationControlInput[5] = getMethodData().getCFL()->getCFLValue();
  
  _lastResidual = currResidual;
  if(_maxResidual < currResidual) {
    _maxResidual = currResidual;
  }
  
  _dissipationControlFunction.evaluate
    (0, _dissipationControlInput, _dissipationControlCoeff);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_FluxSplitter::computeLeftJacobian()
{
  if (!isNull()) {
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    if (m_fluxDerivative.size() != nbEqs) {  
      throw Common::NotImplementedException 
	(FromHere(),"FVMCC_FluxSplitter::computeLeftJacobian()");
    }
    
    _tmpJacobMatrix = 0.;
    for (CFuint iVar = 0; iVar < m_fluxDerivative.size(); ++iVar) {
      // compute d(F)/d(p) where p is each variable of the update state vector 
      // this corresponds to a column in the jacobian matrix 
      (this->*m_fluxDerivative[iVar])(LEFT, iVar, &_tmpJacobMatrix[iVar*nbEqs]);
    }
    _tmpJacobMatrix.transpose(_lFluxJacobian);
    CFLog(DEBUG_MIN, "FVMCC_FluxSplitter::computeLeftJacobian() => _lFluxJacobian = " <<_tmpJacobMatrix << "\n");
  }
}
 
//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_FluxSplitter::computeRightJacobian() 
{
  if (!isNull()) {
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
  if (m_fluxDerivative.size() != nbEqs) {  
    throw Common::NotImplementedException 
      (FromHere(),"FVMCC_FluxSplitter::computeRightJacobian()");
  }
  
  _tmpJacobMatrix = 0.;
  for (CFuint iVar = 0; iVar < m_fluxDerivative.size(); ++iVar) {
    // compute d(F)/d(p) where p is each variable of the update state vector 
    (this->*m_fluxDerivative[iVar])(RIGHT, iVar, &_tmpJacobMatrix[iVar*nbEqs]);
  } 
  _tmpJacobMatrix.transpose(_rFluxJacobian);
  CFLog(DEBUG_MIN, "FVMCC_FluxSplitter::computeLeftJacobian() => _rFluxJacobian = " <<_tmpJacobMatrix << "\n");
  }
}
 
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
