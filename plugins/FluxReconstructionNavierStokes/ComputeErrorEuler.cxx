#include "Framework/MethodCommandProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/ComputeErrorEuler.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider<
    ComputeErrorEuler,FluxReconstructionSolverData,FluxReconstructionNavierStokesModule >
  ComputeErrorEulerProvider("ComputeErrorEuler");

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorEuler::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorEuler::ComputeErrorEuler(const std::string& name) :
  FluxReconstructionSolverCom(name),
  m_eulerVarSet(CFNULL),
  socket_states("states")
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorEuler::~ComputeErrorEuler()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  ComputeErrorEuler::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorEuler::execute()
{
  // get some data from the physical model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal Tt = 0.00365795/m_eulerVarSet->getModel()->getTempRef();
  const CFreal pt = 1.186212306/m_eulerVarSet->getModel()->getPressRef();
  const CFreal M = 0.5;
  const CFreal p = pt/pow(1.0+(gamma-1.0)/2.0*M*M, gamma/(gamma-1.0));
  const CFreal T = Tt/(1.0+(gamma-1.0)/2.0*M*M);
  const CFreal rho = p/(R*T);
  const CFreal a = pow(gamma*p/rho,0.5);
  const CFreal s = p*pow(rho,-gamma);
  const CFreal ht = gamma/(gamma-1.0)*R*T+a*a*M*M/2.0;
  CFreal es = 0.0;
  CFreal eht = 0.0;
  
  // get datahandle for states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  const CFuint nbStates = states.size();

  // loop over the states
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    // dereference states
    State& state   = *states[iState];

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(state,m_solPhysData);

    const CFreal M2 = m_solPhysData[EulerTerm::V]/m_solPhysData[EulerTerm::A];
    const CFreal p2 = m_solPhysData[EulerTerm::P];
    const CFreal T2 = m_solPhysData[EulerTerm::T];
    const CFreal rho2 = m_solPhysData[EulerTerm::RHO];
    const CFreal V2 = m_solPhysData[EulerTerm::V];
    const CFreal s2 = p2*pow(rho2,-gamma);
    const CFreal ht2 = gamma/(gamma-1.0)*R*T2+V2*V2/2.0;
    
    const CFreal errorht = ht - ht2;
    const CFreal errors = s - s2;
    eht += fabs(errorht);
    es += fabs(errors);


  }
  const CFreal errorEnthalpy = eht/(nbStates);
  const CFreal errorEntropy = es/(nbStates);
  CFLog(NOTICE, "error ht: " << errorEnthalpy/ht << "\n");
  CFLog(NOTICE, "error s: " << errorEntropy/s << "\n");
  CFLog(NOTICE, "ht: " << ht << "\n");
  CFLog(NOTICE, "s: " << s << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorEuler::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  FluxReconstructionSolverCom::setup();

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in ComputeErrorEuler!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

