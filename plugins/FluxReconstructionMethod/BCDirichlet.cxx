#include "Framework/MethodStrategyProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BCDirichlet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCDirichlet,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionModule >
  BCDirichletProvider("Dirichlet");

//////////////////////////////////////////////////////////////////////////////

void BCDirichlet::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

BCDirichlet::BCDirichlet(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_inputState()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_functions = vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);

  m_inputVarStr = "";
  setParameter("InputVar",&m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

BCDirichlet::~BCDirichlet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichlet::computeGhostStates(const vector< State* >& intStates,
                                     vector< State* >& ghostStates,
                                     const std::vector< RealVector >& normals,
                                     const std::vector< RealVector >& coords)
{

// Current time
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();  
  CFreal time = subSysStatus->getCurrentTimeDim();

  // number of states
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  RealVector space_time(coords[0].size()+1);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intSol   = *intStates  [iState];
    State& ghostSol = *ghostStates[iState];

 // compute input variables depending on space and time
      for (CFuint i=0; i<coords[iState].size(); ++i)
  	space_time[i]=coords[iState][i];  
      space_time[coords[iState].size()]= time; // time
      m_vFunction.evaluate(space_time, *m_inputState);

  
    // transform to update variables
    State dimState = *m_inputToUpdateVar->transform(m_inputState);

    // adimensionalize the variables if needed and store
    m_varSet->setAdimensionalValues(dimState,ghostSol);

    // modify to ensure that the required value is obtained in the flux point
    ghostSol *= 2.;
    ghostSol -= intSol;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichlet::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                        std::vector< std::vector< RealVector* > >& ghostGrads,
                                        const std::vector< RealVector >& normals,
                                        const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichlet::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // flux point coordinates required
  m_needsSpatCoord = true;

  // set maxNbStatesInCell for variable transformer
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);

  // create state for input from function
  m_inputState = new State();

  // get updateVarSet
  m_varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichlet::unsetup()
{
  CFAUTOTRACE;

  // create state for input from function
  deletePtr(m_inputState);

  // unsetup of the parent class
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichlet::configure ( Config::ConfigArgs& args )
{
  BCStateComputer::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  // get the name of the update variable set
  std::string _updateVarStr = getMethodData().getUpdateVarStr();

  // create the transformer from input to update variables
  if (m_inputVarStr.empty()) {
    m_inputVarStr = _updateVarStr;
  }

  std::string provider =
    VarSetTransformer::getProviderName(physModel->getNameImplementor(), m_inputVarStr, _updateVarStr);

  m_inputToUpdateVar =
    Environment::Factory<VarSetTransformer>::getInstance()
    .getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_inputToUpdateVar.isNotNull());

  // parsing the functions that the user inputed
  m_vFunction.setFunctions(m_functions);

  // brutal fix not to have to change all testcases
  if (m_vars.size() == 2 && m_vars[0] == "x" && m_vars[1] == "y") {m_vars.push_back("t");}
  if (m_vars.size() == 3 && m_vars[0] == "x" && m_vars[1] == "y" && m_vars[2] == "z") {m_vars.push_back("t");}
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
