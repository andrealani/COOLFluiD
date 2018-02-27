#include <numeric>

#include "FluctSplit/FunctionSourceTerm.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FunctionSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitModule>
functionSTProvider("FunctionST");

//////////////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}
      
//////////////////////////////////////////////////////////////////////////////

FunctionSourceTerm::FunctionSourceTerm(const std::string& name) :
  ComputeSourceTermFSM(name),
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
      
//////////////////////////////////////////////////////////////////////////////

FunctionSourceTerm::~FunctionSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::setup()
{
  ComputeSourceTermFSM::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _input.resize(nbEqs);
  _inputVars.resize(nbEqs + dim);
}

//////////////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
					  RealVector& source,
					  const FluctSplit::InwardNormalsData& normalsData)
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const vector<State*>& states = *cell->getStates();
  const CFuint nbStates = states.size();

  _inputVars = 0.;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const State& currState = *(states[iState]);
    const RealVector& coord = currState.getCoordinates();
    // local coordinates
    for (CFuint i = 0; i < dim; ++i) {
      _inputVars[i] += coord[i];
    }
    // local state vector
    for (CFuint i = 0; i < nbEqs; ++i) {
      _inputVars[i+dim] += currState[i];
    }
  }
  _inputVars /= (CFreal)nbStates;
  
  _vFunction.evaluate(_inputVars, _input);
  
  const CFreal volume = socket_volumes.getDataHandle()[cell->getID()];
  source = _input*volume;
  
  CFLog(DEBUG_MAX, "source[" << cell->getID() << "] = " << source << ", vol = " << volume << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
					  std::vector<RealVector>& source,
					  const FluctSplit::InwardNormalsData& normalsData)
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const vector<State*>& states = *cell->getStates();
  const CFuint nbStates = states.size();
  const CFreal volume = socket_volumes.getDataHandle()[cell->getID()];
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const State& currState = *(states[iState]);
    const RealVector& coord = currState.getCoordinates();
    // local coordinates
    for (CFuint i = 0; i < dim; ++i) {
      _inputVars[i] = coord[i];
    }
    // local state vector
    for (CFuint i = 0; i < nbEqs; ++i) {
      _inputVars[i+dim] = currState[i];
    }
    
    _vFunction.evaluate(_inputVars, _input);
    
    source[iState] = _input*volume; // AL: what about this volume? is it correct to use it here?
    
    CFLog(DEBUG_MAX, "source[iState = " << iState << "] = " << source[iState] << ", vol = " << volume << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void FunctionSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ComputeSourceTermFSM::configure(args);
  
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

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
