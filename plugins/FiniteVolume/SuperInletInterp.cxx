#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/SuperInletInterp.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/OldLookupTable.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInletInterp,
                      CellCenterFVMData,
		      FiniteVolumeModule>
superInletInterpFVMCCProvider("SuperInletInterp");

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterp::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< string >
     ("InputFileName","Input data file from which all variables will be extrapolated");
   options.addConfigOption< string >("InputVar","Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletInterp::SuperInletInterp(const std::string& name) :
  FVMCC_BC(name),
  m_inputToUpdateVar(),
  m_lookupState(),
  m_tstate(CFNULL), 
  m_bstate(CFNULL)
{
  addConfigOptionsTo(this);
  
  m_infile = "";
  setParameter("InputFileName",&m_infile); 
  
  m_inputVarStr = "Null";
  setParameter("InputVar",&m_inputVarStr);
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletInterp::~SuperInletInterp()
{
  deletePtr(m_tstate);
  deletePtr(m_bstate);
  
  for (CFuint i = 0; i < m_lookupState.size(); ++i) {
    deletePtr(m_lookupState[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterp::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
  
  const CFreal yCoord = 0.5*(ghostState.getCoordinates()[YY] +
			     innerState.getCoordinates()[YY]);
  
  const CFuint nbEqs = innerState.size();
  for (CFuint i = 0; i < nbEqs; ++i) {
    // interpolated state value in input variables
    (*m_tstate)[i] = m_lookupState[i]->get(yCoord);
  }
  *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
  
  // AL: gory fix just to test
  // if (yCoord > 0.075) (*m_bstate)[11] = std::max((*m_bstate)[11], 50.);
  // (*m_bstate)[11] = std::max((*m_bstate)[11], (CFreal)50.);  
  
  ghostState = 2.*(*m_bstate) - innerState;
  
  for (CFuint i = 0; i < ghostState.size(); ++i) {
    // AL: this gory fix meant for air-11 must be eliminated !!!
    if ((i < 11 || i > 12) && ghostState[i] < 0.) {
      ghostState[i] = (*m_bstate)[i];
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void SuperInletInterp::setup()
{
  CFLog(VERBOSE, "SuperInletInterp::setup() => start\n");
  
  FVMCC_BC::setup();

  CFLog(VERBOSE, "SuperInletInterp::setup() => start1\n");
  
  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  // create the transformer from input to update variables
  if (m_inputVarStr == "Null") {
    m_inputVarStr = getMethodData().getUpdateVarStr();
  }
  
    
  std::string provider = "Identity";
  if (m_inputVarStr != getMethodData().getUpdateVarStr()) {
    cout << m_inputVarStr << " != " << getMethodData().getUpdateVarStr() << endl;
    provider = VarSetTransformer::getProviderName
      (physModel->getConvectiveName(), m_inputVarStr, getMethodData().getUpdateVarStr());
  }
  
  CFLog(VERBOSE, "SuperInletInterp::setup() => start2 << " << provider << "\n");
  
  m_inputToUpdateVar = 
    Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->
    create(physModel->getImplementor());
   
  CFLog(VERBOSE, "SuperInletInterp::setup() => start2 after \n");
  cf_assert(m_inputToUpdateVar.isNotNull());
  
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);
  
  CFLog(VERBOSE, "SuperInletInterp::setup() => start3\n");
  
  m_tstate = new State();
  m_bstate = new State();
  
  if (m_infile != "") {
    if (PE::GetPE().IsParallel()) {
      PE::GetPE().setBarrier();
      for (CFuint p = 0; p < PE::GetPE().GetProcessorCount(); ++p) {
	if (p == PE::GetPE().GetRank ()) {
	  fillTable();
	}
	PE::GetPE().setBarrier();
      }
    }
    else {
      fillTable();
    }
  }
  else {
    CFLog(WARN, "WARNING: SuperInletInterp::setup() => filename not specified!\n");
  } 
  
  CFLog(VERBOSE, "SuperInletInterp::setup() => start4\n");
  
  CFLog(VERBOSE, "SuperInletInterp::setup() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterp::fillTable()
{
  boost::filesystem::path filepath = Environment::DirPaths::getInstance().
    getWorkingDir() / m_infile;
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(filepath);
  
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

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
