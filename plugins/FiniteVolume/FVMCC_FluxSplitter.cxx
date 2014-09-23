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
  _faceArea(0.),
  _firstResidual(0.0),
  _lastResidual(0.0),
  _maxResidual(0.0),
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
  
  if (_dissipationControlDef != "") {
    cf_assert(_dissipationControlVars.size() > 0);
    // i,r,ri,rl,rmax,cfl
    _dissipationControlInput.resize(_dissipationControlVars.size());
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_FluxSplitter::computeFlux(RealVector& result)
{
  const CFuint faceID = getMethodData().getCurrentFace()->getID();
  _faceArea = socket_faceAreas.getDataHandle()[faceID];
  
  // integrate the flux computed at the quadrature points
  (!getMethodData().useAnalyticalConvJacob()) ? 
    integrateFluxOnly(result) : integrateFluxAndJacob(result);
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

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
