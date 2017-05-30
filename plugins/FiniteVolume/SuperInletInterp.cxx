#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/SuperInletInterp.hh"
#include "Framework/StateInterpolator.hh"
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
  options.addConfigOption< string >("InputVar","Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletInterp::SuperInletInterp(const std::string& name) :
  FVMCC_BC(name),
  m_inputToUpdateVar(),
  m_tstate(CFNULL), 
  m_bstate(CFNULL)
{
  addConfigOptionsTo(this);
  
  m_inputVarStr = "Null";
  setParameter("InputVar",&m_inputVarStr);
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletInterp::~SuperInletInterp()
{
  deletePtr(m_tstate);
  deletePtr(m_bstate);
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterp::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
  
  const CFreal yCoord = 0.5*(ghostState.getCoordinates()[YY] +
			     innerState.getCoordinates()[YY]);
  
  SafePtr<StateInterpolator> interp = getStateInterpolator();
  
  const CFuint nbEqs = innerState.size();
  for (CFuint i = 0; i < nbEqs; ++i) {
    // interpolated state value in input variables
    interp->interpolate(i, yCoord, (*m_tstate)[i]);
  }
  *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
  
  // AL: gory fix just to test
  // if (yCoord > 0.075) (*m_bstate)[11] = std::max((*m_bstate)[11], 50.);
  // (*m_bstate)[11] = std::max((*m_bstate)[11], (CFreal)50.);  
  
  ghostState = 2.*(*m_bstate) - innerState;
  
  // for (CFuint i = 0; i < ghostState.size(); ++i) {
  //   // AL: this gory fix meant for air-11 must be eliminated !!!
  //   if ((i < 11 || i > 12) && ghostState[i] < 0.) {
  //     ghostState[i] = (*m_bstate)[i];
  //   }
  // }
}
      
//////////////////////////////////////////////////////////////////////////////

void SuperInletInterp::setup()
{
  CFLog(VERBOSE, "SuperInletInterp::setup() => start\n");
  
  FVMCC_BC::setup();

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  // create the transformer from input to update variables
  if (m_inputVarStr == "Null") {
    m_inputVarStr = getMethodData().getUpdateVarStr();
  }
  
  std::string provider = "Identity";
  if (m_inputVarStr != getMethodData().getUpdateVarStr()) {
    CFLog(VERBOSE, "SuperInletInterp::setup() => inputVar (" << m_inputVarStr 
	  << ") != updateVar (" << getMethodData().getUpdateVarStr() << ")\n");
    provider = VarSetTransformer::getProviderName
      (physModel->getConvectiveName(), m_inputVarStr, getMethodData().getUpdateVarStr());
  }
  
  m_inputToUpdateVar = 
    FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, provider)->
    create(physModel->getImplementor());
  
  cf_assert(m_inputToUpdateVar.isNotNull());
  
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);
  
  m_tstate = new State();
  m_bstate = new State();
  
  CFLog(VERBOSE, "SuperInletInterp::setup() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
