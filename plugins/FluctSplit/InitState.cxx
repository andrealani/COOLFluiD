#include "Common/CFLog.hh"
#include "Common/CFPrintContainer.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/InitState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitState, FluctuationSplitData, FluctSplitModule>
initStateProvider("InitState");

//////////////////////////////////////////////////////////////////////////////

void InitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
  options.addConfigOption< bool >("AdimensionalValues","Input adimensional values.");
  options.addConfigOption< std::vector<CFuint> >
    ("InteractiveVarIDs", "IDs of the variables that will be changed interactively.");
  options.addConfigOption< CFreal, DynamicOption< ValidateOption < PositiveLessThanOne > > >
    ("InteractiveFactor", "Factor to multiply the selected InteractiveVarIDs (should be < 1).");
  options.addConfigOption< std::vector<std::string> >
    ("InitVars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >
    ("InitDef","Definition of the Functions.");
  options.addConfigOption< bool >("UseBlasiusInflow","Generates a Blasius inflow profile based on the distance to the wall, Reynolds number and upstream position");
  options.addConfigOption< CFreal >("ReferenceLength","Reference Length of the system [m].");
  options.addConfigOption< CFreal >("UpstreamPosition","Position [m] aft of the leading edge at which the Blasius profile should be inscribed.");
  options.addConfigOption< CFreal >("ReynoldsNumber","Reynolds number based on reference length, freestream velocity and viscosity.");
  options.addConfigOption< CFreal >("MachNumber","Mach number based on reference length, freestream velocity and viscosity.");
  options.addConfigOption< CFreal >("Gamma","Adiabatic exponent");
}

//////////////////////////////////////////////////////////////////////////////

InitState::InitState(const std::string& name) :
  FluctuationSplitCom(name),
  socket_states("states"),
  _varSet(CFNULL),
  _inputToUpdateVar(),
  _input(CFNULL),
  _tmpFun(),
  _tmpVars(),
  _vFunction(),
  _vInitFunction()
{
  addConfigOptionsTo(this);
   
  _inputAdimensionalValues = false;
  setParameter("AdimensionalValues",&_inputAdimensionalValues);
   
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);
  
  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
  
  _initFunctions = std::vector<std::string>();
  setParameter("InitDef",&_initFunctions);

  _initVars = std::vector<std::string>();
  setParameter("InitVars",&_initVars);
  
  _inputVarStr = "Null";
  setParameter("InputVar",&_inputVarStr);
  
  _interVarIDs = std::vector<CFuint>();
  setParameter("InteractiveVarIDs",&_interVarIDs);
  
  _interFactor = 1.0;
  setParameter("InteractiveFactor",&_interFactor);
  
  m_useBlasius = false;
  setParameter("UseBlasiusInflow",&m_useBlasius);
  
  m_ReferenceLength = 1.0;
  setParameter("ReferenceLength",&m_ReferenceLength);
  
  m_UpstreamPosition = 1.0;
  setParameter("UpstreamPosition",&m_UpstreamPosition);
  
  m_ReynoldsNumber = 1.0;
  setParameter("ReynoldsNumber",&m_ReynoldsNumber);
  
  m_MachNumber = 1.0;
  setParameter("MachNumber",&m_MachNumber);
  
  m_gamma = 1.4;
  setParameter("Gamma",&m_gamma);
}

//////////////////////////////////////////////////////////////////////////////

InitState::~InitState()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitState::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "InitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  Common::SafePtr<std::vector<CFuint> > statesIdx = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbAddVars = _initFunctions.size();
  cf_assert(_tmpVars.size() == dim + nbAddVars);
    
  State dimState;
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    State *const state = states[(*statesIdx)[iState]];
    RealVector& coord = state->getCoordinates();
    if (nbAddVars > 0) {
      // first evaluate expressions for additional variables f1, f2, ...
      // as functions of (x,y,z)
      _vInitFunction.evaluate(coord, _tmpFun);
    }
    
    // set the first dim components of the input vars for the vFunction
    _tmpVars.slice(0, dim) = coord.slice(0, dim);
    
    // set the following additional variables with the result of the
    // already evaluated expressions
    if (nbAddVars > 0) {
      _tmpVars.slice(dim, nbAddVars) = _tmpFun.slice(0, nbAddVars);
    }
    
    // evaluate the state variables as functions of (x,y,z, f1, f2, ...)
    _vFunction.evaluate(_tmpVars, *_input);
    
    // if some interactive variable IDs are specified, multiply
    // those variables by the given factor
    if (_interVarIDs.size() > 0) {
      for (CFuint i = 0; i < _interVarIDs.size(); ++i) {
	(*_input)[_interVarIDs[i]] *= _interFactor;
      }
    }
    
    if(_inputAdimensionalValues) {
      /// @TODO gory fix for now
      _inputToUpdateVar->setLocalID(state->getLocalID());
      *state = *_inputToUpdateVar->transform(_input);
    }
    else {
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *state);
    }
    
   //  if(m_useBlasius)
//     {
//       if(coord[0]<0.0 || coord[1]<0.0) std::cout<<"ERROR negative gridpoint ("<<coord[0]<<","<<coord[1]<<","<<coord[2]<<")"<<std::endl;
//       CFreal y1		=	0.0;
//       CFreal y2		=	0.0;
//       CFreal y3		=	0.4696;
//       CFreal eta 	=       0.0;
//       CFreal h		=	0.01;	
//       CFreal eta_end	=	coord[1]/sqrt(2.0/m_ReynoldsNumber*(m_UpstreamPosition/m_ReferenceLength+coord[0]));
//       CFreal y1temp,y2temp,y3temp;
//       CFreal k11,k12,k13;
//       CFreal k21,k22,k23;
//       CFreal k31,k32,k33;
//       CFreal k41,k42,k43;
      
//       if(eta_end>10.0) y2=1.00;
//       else
//       {
// 	while(eta<eta_end)
// 	{
// 	  k11=h*y2; k12=h*y3; k13	=h*(-y3*y1);
// 	  y1temp = y1+0.5*k11; y2temp = y2+0.5*k12; y3temp= y3+0.5*k13;
	  
// 	  k21=h*y2temp; k22=h*y3temp; k23 = h*(-y1temp*y3temp);
// 	  y1temp= y1+0.5*k21; y2temp = y2+0.5*k22; y3temp = y3+0.5*k23;
	  
// 	  k31=h*y2temp; k32=h*y3temp; k33=h*(-y1temp*y3temp);
// 	  y1temp= y1+k31; y2temp = y2+k32; y3temp = y3+k33;
	  
// 	  k41=h*y2temp; k42=h*y3temp; k43=h*(-y1temp*y3temp);
	  
// 	  y1=y1+(k11+2.0*k21+2.0*k31+k41)/6.0;
// 	  y2=y2+(k12+2.0*k22+2.0*k32+k42)/6.0;
// 	  y3=y3+(k13+2.0*k23+2.0*k33+k43)/6.0;
// 	  eta = eta + h;
// 	}
//       }
//       blasius_inflow=y2;
//       blasius_energy_inflow=1.0/m_gamma/(m_gamma-1.0)/m_MachNumber/m_MachNumber;
//       ///std::cout<<"y="<<coord[1]<<" eta="<<eta-h<<" f'(eta)="<<y2<<" E-k="<<blasius_energy_inflow[state]<<std::endl;
//       (*state)[1]=(*state)[1]*blasius_inflow;
//       (*state)[4]=(*state)[0]*blasius_energy_inflow+0.5/(*state)[0]*(*state)[1]*(*state)[1];
//       (*state)[2]=0.0;
//       (*state)[3]=0.0;
//     }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
InitState::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void InitState::setup()
{
  CFAUTOTRACE;
  
  FluctuationSplitCom::setup();

  _input = new State();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  if (_initFunctions.size() > 0) {
    _tmpFun.resize(_initFunctions.size());
    cf_assert(_tmpFun.size() > 0);
  }
  
  _tmpVars.resize(dim + _tmpFun.size());
  cf_assert(_tmpVars.size() > 0);
  
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  
  _inputToUpdateVar->setup(maxNbStatesInCell);
  _varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::unsetup()
{
  CFAUTOTRACE;

  deletePtr(_input);

  FluctuationSplitCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FluctuationSplitCom::configure(args);

 std::string name = getMethodData().getNamespace();
 Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
   (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  // create the transformer from input to update variables
  if (_inputVarStr == "Null") {
    _inputVarStr = getMethodData().getUpdateVarStr();
  }
  
 std::string provider = VarSetTransformer::getProviderName
    (physModel->getNameImplementor(), _inputVarStr, getMethodData().getUpdateVarStr());
  
  // fix to allow identity transformations
  if ( ! Environment::Factory<VarSetTransformer>::getInstance().exists (provider) )
    provider = "Identity";
  
  _inputToUpdateVar =
    Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->create(physModel->getImplementor());
  
 // cout << CFPrintContainer<vector<string> >("_functions  = ", &_functions) << endl;
 // cout << CFPrintContainer<vector<string> >("_vars  = ", &_vars) << endl;
  
  cf_assert(_inputToUpdateVar.isNotNull());
  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // rethrow the exception to signal the error to the user
  } 
  
 // cout << CFPrintContainer<vector<string> >("_initFunctions  = ", &_initFunctions) << endl;
 // cout << CFPrintContainer<vector<string> >("_initVars  = ", &_initVars) << endl;
  
  if (_initFunctions.size() > 0) {
    _vInitFunction.setFunctions(_initFunctions);
    _vInitFunction.setVariables(_initVars);
    try {
      _vInitFunction.parse();
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // rethrow the exception to signal the error to the user
    }  
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
