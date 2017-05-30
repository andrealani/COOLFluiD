#include "FiniteVolume/FiniteVolume.hh"
#include "SpecialSuperInlet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SpecialSuperInlet, CellCenterFVMData, FiniteVolumeModule> SpecialSuperInletFVMCCProvider("SpecialSuperInletFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SpecialSuperInlet::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
  options.addConfigOption< bool >("AdimensionalValues","Flag to input adimensional values.");
}

//////////////////////////////////////////////////////////////////////////////

SpecialSuperInlet::SpecialSuperInlet(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _bCoord(),
  _dimState(CFNULL),
  _inputToUpdateVar(),
  _input(CFNULL)
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

   _inputVarStr = "Null";
   setParameter("InputVar",&_inputVarStr);

  _inputAdimensionalValues = false;
   setParameter("AdimensionalValues",&_inputAdimensionalValues);

}

//////////////////////////////////////////////////////////////////////////////

SpecialSuperInlet::~SpecialSuperInlet()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpecialSuperInlet::setup()
{

  FVMCC_BC::setup();

  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 2);

  _varSet = getMethodData().getUpdateVar();

  _dimState = new State();

  _input = new State();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _inputToUpdateVar->setup(maxNbStatesInCell);

}

//////////////////////////////////////////////////////////////////////////////

void SpecialSuperInlet::unsetup()
{

  deletePtr(_input);
  deletePtr(_dimState);

  FVMCC_BC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SpecialSuperInlet::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);

  const std::string name = getMethodData().getNamespace();

  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);

  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  /// @todo Remove this !!!
  _updateVarStr = getMethodData().getUpdateVarStr();

  // create the transformer from input to update variables
  if (_inputVarStr == "Null") {
    _inputVarStr = _updateVarStr;
  }

  const std::string provider = VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), _inputVarStr, _updateVarStr);
// CFout << "Trying to use provider: " << provider <<"\n";
  _inputToUpdateVar =
      FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, provider)->
    create(physModel->getImplementor());
  
  cf_assert(_inputToUpdateVar.isNotNull());

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

void SpecialSuperInlet::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  for(CFuint iDim=0; iDim < nbDim; ++iDim)
  {
    _variables[iDim] = _bCoord[iDim];
  }
  _variables[nbDim] = SubSystemStatusStack::getActive()->getNbIter();
  _variables[nbDim+1] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  //Evaluate the function
  _vFunction.evaluate(_variables,*_input);

// CFout << "Input: " << *_input <<"\n";
  //Set the state value
  if (_inputAdimensionalValues){
    *ghostState = *_inputToUpdateVar->transform(_input);
  }
  else{
    *_dimState = *_inputToUpdateVar->transform(_input);
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
  }
  
  *ghostState *= 2.;
  *ghostState -= *innerState;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
