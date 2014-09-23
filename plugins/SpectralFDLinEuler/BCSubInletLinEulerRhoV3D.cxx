#include "Framework/MethodStrategyProvider.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "LinEuler/LinEuler3DVarSet.hh"
#include "LinEuler/LinEulerTerm.hh"

#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"
#include "SpectralFDLinEuler/BCSubInletLinEulerRhoV3D.hh"

#include "Common/NotImplementedException.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::LinearizedEuler;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSubInletLinEulerRhoV3D,SpectralFDMethodData,BCStateComputer,SpectralFDLinEulerModule >
    BCSubInletLinEulerRhoV3DProvider("SubInletLinEulerRhoV3D");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BCSubInletLinEulerRhoV3D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_meanflow);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletLinEulerRhoV3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletLinEulerRhoV3D::BCSubInletLinEulerRhoV3D(const std::string& name) :
  BCStateComputer(name),
  socket_meanflow("meanflow"),
  m_linEulerVarSet(CFNULL),
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

BCSubInletLinEulerRhoV3D::~BCSubInletLinEulerRhoV3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletLinEulerRhoV3D::computeGhostStates(const vector< State* >& intStates,
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

  cf_assert(nbrStates <= m_extraVars->size());

  // create real vector for space and time independent variables
  RealVector space_time(coords[0].size()+1); 

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    // check dimensions of internal and ghost states
    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);
    cf_assert((*m_extraVars)[iState]->size() == 5);

     /// acquaintance of the PhysicalModel
     Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model =  m_linEulerVarSet->getModel();

     CFuint stateID = intState.getLocalID();

     RealVector meanflow_state = _model->getMeanFlowState(stateID);

    // set the extra variables
    m_linEulerVarSet->setExtraPhysicalVars(&meanflow_state);


    // compute input variables depending on space and time

    // const CFreal dt = SubSystemStatusStack::getActive()->getDT();

    for (CFuint i=0; i<coords[iState].size(); ++i)
        space_time[i]=coords[iState][i];
                                                             //std::cout<< dt << std::endl;
        space_time[coords[iState].size()]= time; // time

        m_vFunction.evaluate(space_time,m_inputVars);
                                                             //std::cout<< m_inputVars << std::endl;

 /*   const CFreal rho0   = m_intSolPhysData[LinEulerTerm::rho0 ];
    const CFreal U0     = m_intSolPhysData[LinEulerTerm::U0   ];
    const CFreal V0     = m_intSolPhysData[LinEulerTerm::V0   ];*/
/*

    // non-dimensionalize input variables
    m_inputVars[0] /= m_rhoRef;
    m_inputVars[1] /= m_velRef;
    m_inputVars[2] /= m_velRef;
    m_inputVars[3] /= m_velRef;    */

    // set the physical data starting from the inner state
    m_linEulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // set ghost state quantities
    const CFreal rhoGhost = 2.0*m_inputVars[0] - m_intSolPhysData[LinEulerTerm::rho];
    const CFreal vxGhost  = 2.0*m_inputVars[1] - m_intSolPhysData[LinEulerTerm::u ];
    const CFreal vyGhost  = 2.0*m_inputVars[2] - m_intSolPhysData[LinEulerTerm::v ];
    const CFreal vzGhost  = 2.0*m_inputVars[3] - m_intSolPhysData[LinEulerTerm::w ];
    const CFreal pGhost   = m_intSolPhysData[LinEulerTerm::p];

    //set the physical data for the ghost state
    m_ghostSolPhysData[LinEulerTerm::p]   = pGhost;
    m_ghostSolPhysData[LinEulerTerm::rho] = rhoGhost;
    m_ghostSolPhysData[LinEulerTerm::u]  = vxGhost;
    m_ghostSolPhysData[LinEulerTerm::v]  = vyGhost;
    m_ghostSolPhysData[LinEulerTerm::w]  = vzGhost;

    // set the ghost state from its physical data
    m_linEulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletLinEulerRhoV3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                               std::vector< std::vector< RealVector* > >& ghostGrads,
                                               const std::vector< RealVector >& normals,
                                               const std::vector< RealVector >& coords)
{
  throw Common::NotImplementedException(FromHere(),"BCSubInletLinEulerRhoV3D::computeGhostGradients(): this should not be needed for linearized Euler!!!");
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletLinEulerRhoV3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // flux point coordinates required
  m_needsSpatCoord = true;

  // extra variables required
  m_needsExtraVars = true;

  // get linearized Euler 3D varset
  m_linEulerVarSet = getMethodData().getUpdateVar().d_castTo<LinEuler3DVarSet>();
  if (m_linEulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException(FromHere(),"Update variable set is not LinEuler3DVarSet in BCSubInletLinEulerRhoV3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_linEulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_linEulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // resize m_inputVars
  m_inputVars.resize(4);


// get reference values for non-dimensonalization
  RealVector refPhysData =  m_linEulerVarSet->getModel()->getReferencePhysicalData();
  m_rhoRef = refPhysData[LinEulerTerm::rho];
  m_velRef = refPhysData[LinEulerTerm::u  ];


  //Get mean flow socket
  socket_meanflow.setParentNamespace( getMethodData().getNamespace() );

}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletLinEulerRhoV3D::configure ( Config::ConfigArgs& args )
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
