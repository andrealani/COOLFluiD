#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/BCSubInletEulerRhoV3D.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSubInletEulerRhoV3D,SpectralFDMethodData,BCStateComputer,SpectralFDNavierStokesModule >
  BCSubInletEulerRhoV3DProvider("SubInletEulerRhoV3D");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerRhoV3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletEulerRhoV3D::BCSubInletEulerRhoV3D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_inputVars(),
  m_rhoRef(),
  m_velRef()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_functions = vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletEulerRhoV3D::~BCSubInletEulerRhoV3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerRhoV3D::computeGhostStates(const vector< State* >& intStates,
                                               vector< State* >& ghostStates,
                                               const std::vector< RealVector >& normals,
                                               const std::vector< RealVector >& coords)
{
  // Current time
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();  
  CFreal time = subSysStatus->getCurrentTimeDim();
  
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some data from the physical model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  RealVector space_time(coords[0].size()+1);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

      // compute input variables depending on space and time
      for (CFuint i=0; i<coords[iState].size(); ++i)
  	space_time[i]=coords[iState][i];  
      space_time[coords[iState].size()]= time; // time
      m_vFunction.evaluate(space_time,m_inputVars);

    // non-dimensionalize input variables
    m_inputVars[0] /= m_rhoRef;
    m_inputVars[1] /= m_velRef;
    m_inputVars[2] /= m_velRef;
    m_inputVars[3] /= m_velRef;

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // set ghost state quantities
    const CFreal rhoGhost = 2.0*m_inputVars[0] - m_intSolPhysData[EulerTerm::RHO];
    const CFreal vxGhost  = 2.0*m_inputVars[1] - m_intSolPhysData[EulerTerm::VX ];
    const CFreal vyGhost  = 2.0*m_inputVars[2] - m_intSolPhysData[EulerTerm::VY ];
    const CFreal vzGhost  = 2.0*m_inputVars[3] - m_intSolPhysData[EulerTerm::VZ ];
    const CFreal pGhost   = m_intSolPhysData[EulerTerm::P];

    //set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::P]   = pGhost;
    m_ghostSolPhysData[EulerTerm::RHO] = rhoGhost;
    m_ghostSolPhysData[EulerTerm::VX]  = vxGhost;
    m_ghostSolPhysData[EulerTerm::VY]  = vyGhost;
    m_ghostSolPhysData[EulerTerm::VZ]  = vzGhost;
    m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*pGhost
                                          + 0.5*rhoGhost*
                                          (vxGhost*vxGhost +
                                           vyGhost*vyGhost +
                                           vzGhost*vzGhost)
                                         )/rhoGhost;

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerRhoV3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

void BCSubInletEulerRhoV3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // flux point coordinates required
  m_needsSpatCoord = true;

  // get Euler 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCSubInletEulerRhoV3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // resize m_inputVars
  m_inputVars.resize(4);

  // get reference values for non-dimensonalization
  RealVector refPhysData =  m_eulerVarSet->getModel()->getReferencePhysicalData();
  m_rhoRef = refPhysData[EulerTerm::RHO];
  m_velRef = refPhysData[EulerTerm::V  ];
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerRhoV3D::configure ( Config::ConfigArgs& args )
{
  BCStateComputer::configure(args);

  // parsing the functions that the user inputed
  m_vFunction.setFunctions(m_functions);
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

  }  // namespace SpectralFD

}  // namespace COOLFluiD
