#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/SuperInlet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInlet, CellCenterFVMData, FiniteVolumeModule> 
superInletFVMCCProvider("SuperInletFVMCC");

MethodCommandProvider<SuperInlet, CellCenterFVMData, FiniteVolumeModule> 
dirichletConditionFVMCCProvider("DirichletConditionFVMCC");
      
//////////////////////////////////////////////////////////////////////////////

void SuperInlet::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
  options.addConfigOption< bool >("AdimensionalValues","Flag to input adimensional values.");
  options.addConfigOption< std::vector<CFuint> >
    ("InteractiveVarIDs", "IDs of the variables that will be changed interactively.");
  options.addConfigOption< CFreal, DynamicOption< ValidateOption < PositiveLessThanOne > > >
  ("InteractiveFactor", "Factor to multiply the selected InteractiveVarIDs (should be < 1).");
}

//////////////////////////////////////////////////////////////////////////////

SuperInlet::SuperInlet(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _bCoord(),
  _xyzIter(),
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
  
  _interVarIDs = std::vector<CFuint>();
  setParameter("InteractiveVarIDs",&_interVarIDs);  
  
  _interFactor = 1.0;
  setParameter("InteractiveFactor",&_interFactor);
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInlet::~SuperInlet()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::setup()
{

  FVMCC_BC::setup();

  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
   
  _xyzIter.resize(PhysicalModelStack::getActive()->getDim() + 1);
  
  _varSet = getMethodData().getUpdateVar();
  
  _dimState = new State();
  
  _input = new State();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _inputToUpdateVar->setup(maxNbStatesInCell);
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::unsetup()
{

  deletePtr(_input);
  deletePtr(_dimState);

  FVMCC_BC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);

  const std::string name = getMethodData().getNamespace();

  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);

  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

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

void SuperInlet::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const bool hasIter = (_vars.size() == dim + 1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;
  
  if (!hasIter) {
    //Evaluate the function
    _vFunction.evaluate(_bCoord,*_input);
  }
  else {
    // [x,y,z] and iteration are fed to the evaluate function
    for (CFuint i = 0; i < dim; ++i) {
      _xyzIter[i] = _bCoord[i]; 
    }
    _xyzIter[dim] = SubSystemStatusStack::getActive()->getNbIter();
    
    //Evaluate the function
    _vFunction.evaluate(_xyzIter,*_input);
  }
  
  // if some interactive variable IDs are specified, multiply 
  // those variables by the given factor
  if (_interVarIDs.size() > 0) {
    for (CFuint i = 0; i < _interVarIDs.size(); ++i) {
      (*_input)[_interVarIDs[i]] *= _interFactor;
    }
  }
  
  // CFout << "Input: " << *_input <<"\n";
  //Set the state value
  if (_inputAdimensionalValues){
    *ghostState = *_inputToUpdateVar->transform(_input);
  }
  else{
    *_dimState = *_inputToUpdateVar->transform(_input);
    
    // CFout << "DimState: " << *_dimState <<"\n";
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
    
    // CFout << "ADimState: " << *ghostState <<"\n";
  }
  *ghostState *= 2.;
  *ghostState -= *innerState;
  
  // CFLog(INFO, "ADimState: " << *ghostState <<"\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
