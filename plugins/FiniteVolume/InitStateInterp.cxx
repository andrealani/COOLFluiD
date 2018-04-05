#include "FiniteVolume/FiniteVolume.hh"

#include "Common/OldLookupTable.hh"
#include "InitStateInterp.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateInterp, CellCenterFVMData, FiniteVolumeModule>
initStateInterpProvider("InitStateInterp");

//////////////////////////////////////////////////////////////////////////////

void InitStateInterp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("InputFileName","Input data file from which all variables will be extrapolated");
  options.addConfigOption< std::string >("InputInterpVar","Input interpolation variables.");
}
      
//////////////////////////////////////////////////////////////////////////////

InitStateInterp::InitStateInterp(const std::string& name) :
  InitState(name),
  m_inputInterpToUpdateVar(),
  m_lookupState(),
  m_tstate(CFNULL), 
  m_bstate(CFNULL)
{
  addConfigOptionsTo(this);
  
  m_infile = "";
  setParameter("InputFileName",&m_infile); 
  
  m_inputInterpVarStr = "Null";
  setParameter("InputInterpVar",&m_inputInterpVarStr);
}

//////////////////////////////////////////////////////////////////////////////

InitStateInterp::~InitStateInterp()
{ 
}

//////////////////////////////////////////////////////////////////////////////

void InitStateInterp::executeOnTrs()
{
  CFLog(VERBOSE, "InitStateInterp::executeOnTrs() => start\n");
  
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitStateInterp::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");
  
  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitStateInterp not applied to InnerFaces!!!");
  }

  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  std::vector<CFuint>::iterator itd;
  if(_inputAdimensionalValues)
  {
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      
      const CFreal yCoord = currState->getCoordinates()[YY];
      const CFuint nbEqs = currState->size();
      for (CFuint i = 0; i < nbEqs; ++i) {
	// interpolated state value in input variables
	(*m_tstate)[i] = m_lookupState[i]->get(yCoord);
      }
      
      *currState = *_inputToUpdateVar->transform(m_tstate);
    }
  }
  else
  {
    State dimState;
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      const CFuint nbEqs = currState->size();
      const CFreal yCoord = currState->getCoordinates()[YY];
      
      //if (currState->getCoordinates()[XX] < -0.05 && yCoord < 0.075) {
      // if (yCoord < 0.075) {
      for (CFuint i = 0; i < nbEqs; ++i) {
	// interpolated state value in input variables
	(*m_tstate)[i] = m_lookupState[i]->get(yCoord);
      }
      CFLog(DEBUG_MIN, "InitStateInterp::executeOnTrs() => m_tstate = " << *m_tstate << "\n");
      dimState = *m_inputInterpToUpdateVar->transform(m_tstate);
      //  }
      //       else {
      // 	_vFunction.evaluate(currState->getCoordinates(), *m_tstate);
      // 	dimState = *_inputToUpdateVar->transform(m_tstate);
      //       }
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
  CFLog(VERBOSE, "InitStateInterp::executeOnTrs() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void InitStateInterp::setup()
{
  CFAUTOTRACE;
  
  InitState::setup();
  
  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  if (m_inputInterpVarStr == "Null") {
    m_inputInterpVarStr = getMethodData().getUpdateVarStr();
  }
  const std::string provider = VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), m_inputInterpVarStr, getMethodData().getUpdateVarStr());
  
  m_inputInterpToUpdateVar =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, provider)->
    create(physModel->getImplementor());
  cf_assert(m_inputInterpToUpdateVar.isNotNull());
  
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputInterpToUpdateVar->setup(maxNbStatesInCell);
  
  m_tstate = new State();
  m_bstate = new State();
  
  if (m_infile != "") {
    if (PE::GetPE().IsParallel()) {
      PE::GetPE().setBarrier(name);
      for (CFuint p = 0; p < PE::GetPE().GetProcessorCount(name); ++p) {
	if (p == PE::GetPE().GetRank(name)) {
	  fillTable();
	}
	PE::GetPE().setBarrier(name);
      }
    }
    else {
      fillTable();
    }
  }
  else {
    cout << "WARNING: InitStateInterp::setup() => filename not specified!"<< endl;
  }
}

//////////////////////////////////////////////////////////////////////////////

void InitStateInterp::fillTable()
{
  boost::filesystem::path filepath = Environment::DirPaths::getInstance().
    getWorkingDir() / m_infile;
  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& fin = (*fhandle)->open(filepath);
  
  string variables;
  // read the first line with the variables names
  getline(fin,variables);
  CFuint nbPoints = 0;
  fin >> nbPoints;
  
  // allocate the look up tables
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_lookupState.resize(nbEqs);
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i] = new LookUpTable<CFreal, CFreal>(nbPoints);
  } 
  
  // nbEqs + "y" AL: hardcoded here
  CFreal ycoord = 0.;
  CFreal tmpVar = 0.;
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
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////
      
void InitStateInterp::unsetup()
{
  CFAUTOTRACE;
  
  deletePtr(m_tstate);
  deletePtr(m_bstate);
  for (CFuint i = 0; i < m_lookupState.size(); ++i) {
    deletePtr(m_lookupState[i]);
  }
  
  InitState::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
