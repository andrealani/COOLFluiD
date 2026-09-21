#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BCDirichletFromFile.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Common/OldLookupTable.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCDirichletFromFile,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionModule >
  BCDirichletFromFileProvider("DirichletFromFile");

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::defineConfigOptions(Config::OptionList& options)
{
   options.template addConfigOption< std::string >("InputFileName", "Input data file from which all variables will be extrapolated.");
   options.template addConfigOption< std::string >("InputVar", "Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

BCDirichletFromFile::BCDirichletFromFile(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_dimState(CFNULL),
  m_input(CFNULL),
  m_lookupState()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);
  
  m_infile = "";
  setParameter("InputFileName",&m_infile); 
  
  m_inputVarStr = "Null";
  setParameter("InputVar",&m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

BCDirichletFromFile::~BCDirichletFromFile()
{
  CFAUTOTRACE;
  
  for (CFuint i = 0; i < m_lookupState.size(); ++i) 
  {
    deletePtr(m_lookupState[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::computeGhostStates(const vector< State* >& intStates,
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
  
  const CFuint nbEqs = intStates[0]->size();

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intSol   = *intStates  [iState];
    State& ghostSol = *ghostStates[iState];
    
    const CFreal yCoord = coords[iState][1];
    for (CFuint i = 0; i < nbEqs; ++i) 
    {
      // interpolated state value in input variables
      (*m_input)[i] = m_lookupState[i]->get(yCoord);
    }
    
    // transform to update variables
    *m_dimState = *m_inputToUpdateVar->transform(m_input);
//CFLog(INFO, "coord: " << yCoord << ", state: " << *m_input << ", transfoState: " << *m_dimState << "\n");
    // adimensionalize the variables if needed and store
    m_varSet->setAdimensionalValues(*m_dimState,ghostSol);

    // modify to ensure that the required value is obtained in the flux point
    ghostSol *= 2.;
    ghostSol -= intSol;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

void BCDirichletFromFile::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                             const std::vector< Framework::State* >& intStates,
                                             const std::vector< Framework::State* >& ghostStates,
                                             const std::vector< RealVector >& unitNormals,
                                             const std::vector< RealVector >& flxPntCoords,
                                             std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrStates = intStates.size();

  if (nbrStates == 0)
  {
    return;
  }

  const CFuint nbrEqs = intStates[0]->size();

  if (m_prescStates.size() < nbrStates)
  {
    m_prescStates.resize(nbrStates);
    m_prescStatePtrs.resize(nbrStates);
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      m_prescStates[iState].resize(nbrEqs);
      m_prescStatePtrs[iState] = &m_prescStates[iState];
    }
  }
  if (m_prescGradVars.nbRows() != nbrEqs || m_prescGradVars.nbCols() < nbrStates)
  {
    m_prescGradVars.resize(nbrEqs,nbrStates);
  }

  // prescribed states at the flux points
  computeBndStates(intStates,ghostStates,unitNormals,flxPntCoords,m_prescStatePtrs);

  // g_b = g(U_prescribed)
  getMethodData().getDiffusiveVar()->setGradientVars(m_prescStatePtrs,m_prescGradVars,nbrStates);

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
    {
      (*bndGradVars[iState])[iEq] = m_prescGradVars(iEq,iState);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::computeBndStates(const std::vector< Framework::State* >& intStates,
                                           const std::vector< Framework::State* >& ghostStates,
                                           const std::vector< RealVector >& unitNormals,
                                           const std::vector< RealVector >& flxPntCoords,
                                           std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();

  // U_b = 0.5*(U + U_ghost), with U_ghost = 2*U_file - U
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    *bndStates[iState] = 0.5*(*intStates[iState] + *ghostStates[iState]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // flux point coordinates required
  m_needsSpatCoord = true;
 
  m_dimState = new State();
  m_input    = new State();
  const CFuint maxNbStatesInCell = Framework::MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_inputToUpdateVar->setup(maxNbStatesInCell);
  
  const std::string nsp = this->getMethodData().getNamespace();
  
  if (m_infile != "") {
    if (Common::PE::GetPE().IsParallel()) {
      Common::PE::GetPE().setBarrier(nsp);
      for (CFuint p = 0; p < Common::PE::GetPE().GetProcessorCount(nsp); ++p) {
	if (p == Common::PE::GetPE().GetRank(nsp)) {
	  fillTable();
	}
	Common::PE::GetPE().setBarrier(nsp);
      }
    }
    else {
      fillTable();
    }
  }
  else {
    CFLog(WARN, "WARNING: SuperInletInterp::setup() => filename not specified!\n");
  }

  // get updateVarSet
  m_varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::unsetup()
{
  CFAUTOTRACE;
  
  deletePtr(m_dimState);
  deletePtr(m_input);

  // unsetup of the parent class
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::configure ( Config::ConfigArgs& args )
{
  BCStateComputer::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  // get the name of the update variable set
  std::string updateVarStr = getMethodData().getUpdateVarStr();

  // create the transformer from input to update variables
  if (m_inputVarStr.empty()) 
  {
    m_inputVarStr = updateVarStr;
  }

  std::string provider = VarSetTransformer::getProviderName(physModel->getNameImplementor(), m_inputVarStr, updateVarStr);

  m_inputToUpdateVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_inputToUpdateVar.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void BCDirichletFromFile::fillTable()
{
  boost::filesystem::path filepath = Environment::DirPaths::getInstance().getWorkingDir() / m_infile;

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  std::ifstream& fin = fhandle->open(filepath);
  
  std::string variables;
  // read the first line with the nb of points
  getline(fin, variables);
  CFLog(VERBOSE, "Variables: " << variables << "\n");
  CFuint nbPoints = 0;
  fin >> nbPoints; 
  CFLog(VERBOSE, "BCDirichletFromFile::fillTable() => nbPoints = " << nbPoints << "\n");
  
  // allocate the look up tables
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  m_lookupState.resize(nbEqs);
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i] = new Common::LookUpTable<CFreal, CFreal>(nbPoints);
  } 
  
  // nbEqs + "y" AL: hardcoded here
  CFreal ycoord;
  CFreal tmpVar;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    fin >> ycoord;
    for (CFuint i = 0; i < nbEqs; ++i) {
      fin >> tmpVar;
      m_lookupState[i]->insert(ycoord, tmpVar);
    }
  }
  
  // sort the data 
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i]->sortKeys();
  }
  
  fin.close();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
