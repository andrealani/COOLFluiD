#include "FiniteVolume/FunctionSourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FunctionSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeModule>
functionSTFVMCCProvider("FunctionST");
      
//////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}
      
//////////////////////////////////////////////////////////////////////////////

FunctionSourceTerm::FunctionSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _input(),
  _inputVars(),
  _vFunction()
{
  addConfigOptionsTo(this);

  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);
  
  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////

FunctionSourceTerm::~FunctionSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _input.resize(nbEqs);
  _inputVars.resize(nbEqs + dim);
}

//////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::unsetup()
{
  CFAUTOTRACE;
  
  ComputeSourceTermFVMCC::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::computeSource(GeometricEntity *const element,
				       RealVector& source,
				       RealMatrix& jacobian)
{
  CFLogDebugMin( "FunctionSourceTerm::computeSource()" << "\n");

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const State& currState = *element->getState(0);
  RealVector& coord = currState.getCoordinates();

  // local coordinates
  for (CFuint i = 0; i < dim; ++i) {
    _inputVars[i] = coord[i];
  }
  // local state vector
  for (CFuint i = 0; i < nbEqs; ++i) {
    _inputVars[i+dim] = currState[i];
  }
  
  _vFunction.evaluate(_inputVars, _input);
  source = _input;

  CFLog(DEBUG_MAX, "FunctionSourceTerm::computeSource() => source = "
	<< source << "\n");
  
  const CFreal volume = socket_volumes.getDataHandle()[element->getID()];
  
  CFLog(DEBUG_MAX, "source[" << element->getID() << "] = " << source << ", vol = " << volume << "\n");
  source *= volume;
}

//////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ComputeSourceTermFVMCC::configure(args);
  
  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
