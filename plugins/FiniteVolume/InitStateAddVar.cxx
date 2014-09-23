

#include "FiniteVolume/FiniteVolume.hh"
#include "InitStateAddVar.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateAddVar, CellCenterFVMData, FiniteVolumeModule>
initStateAddVarProvider("InitStateAddVar");

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >
    ("InitVars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >
    ("InitDef","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

InitStateAddVar::InitStateAddVar(const std::string& name) :
  InitState(name),
  _tmpFun(),
  _tmpVars()
{
  addConfigOptionsTo(this);
  _initFunctions = std::vector<std::string>();
  setParameter("InitDef",&_initFunctions);

  _initVars = std::vector<std::string>();
  setParameter("InitVars",&_initVars);
}

//////////////////////////////////////////////////////////////////////////////

InitStateAddVar::~InitStateAddVar()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::configure ( Config::ConfigArgs& args )
{
  InitState::configure(args);

  _vInitFunction.setFunctions(_initFunctions);
  _vInitFunction.setVariables(_initVars);
  try {
    _vInitFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::setup()
{
  InitState::setup();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _tmpFun.resize(_initFunctions.size());
  cf_assert(_tmpFun.size() > 0);

  _tmpVars.resize(dim + _tmpFun.size());
  cf_assert(_tmpVars.size() > 0);
}

//////////////////////////////////////////////////////////////////////////////

void InitStateAddVar::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitStateAddVar::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");

  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitStateAddVar not applied to InnerFaces!!!");
  }

  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbAddVars = _initFunctions.size();
  cf_assert(_tmpVars.size() == dim + nbAddVars);

  State dimState;
  std::vector<CFuint>::iterator itd;
  for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
    State* const currState = states[(*itd)];
    RealVector& coord = currState->getCoordinates();

    // first evaluate expressions for additional variables f1, f2, ...
    // as functions of (x,y,z)
    _vInitFunction.evaluate(coord, _tmpFun);

    // set the first dim components of the input vars for the vFunction
    _tmpVars.slice(0,dim) = coord.slice(0,dim);

    // set the following additional variables with the result of the
    // already evaluated expressions
    _tmpVars.slice(dim, nbAddVars) = _tmpFun.slice(0,nbAddVars);

    // evaluate the state variables as functions of (x,y,z, f1, f2, ...)
    _vFunction.evaluate(_tmpVars, *_input);

    if(_inputAdimensionalValues)
    {
      *currState = *_inputToUpdateVar->transform(_input);
    }
    else
    {
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
