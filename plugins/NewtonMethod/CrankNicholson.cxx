// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CrankNicholson.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CrankNicholson,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
crankNichsolsonMethodProvider("CrankNicholson");

//////////////////////////////////////////////////////////////////////////////

CrankNicholson::CrankNicholson(const std::string& name)
  : NewtonIterator(name)
{
  // new default commands are set for the Crank-Nicholson scheme
  m_setupStr   = "CrankNichSetup";
  m_unSetupStr = "CrankNichUnSetup";
  m_intermediateStr = "CrankNichIntermediate";
}

//////////////////////////////////////////////////////////////////////////////

CrankNicholson::~CrankNicholson()
{
}

//////////////////////////////////////////////////////////////////////////////

void CrankNicholson::takeStepImpl()
{
  // update subsystemstatus if first subiteration step
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  if(subSysStatus->isSubIterationFirstStep()){
    subSysStatus->updateNbIter();
    subSysStatus->updateTimeStep();
    //    getConvergenceMethodData()->getCFL()->update();
  }
  
  // do a prepare step, usually backing up the solution to
  // pastStates
  CFLogDebugMed( "CrankNicholson::takeStep(): calling Prepare step\n");
  m_prepare->execute();

  // at the first time step, compute and back up spatial rhs of initial solution
  if(subSysStatus->getNbIter() == 1)
  {
    // compute spatial rhs
    getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();

    // here compute only the RHS, not the jacobian matrix
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(false);
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(0.5);
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(true);

    // put the computed spatial rhs in pastRhs
    m_data->setAchieved(true);
    const bool isLastStep = subSysStatus->isSubIterationLastStep();
    subSysStatus->setIsSubIterationLastStep(true);
    m_intermediate->execute();
    subSysStatus->setIsSubIterationLastStep(isLastStep);
  }

  // update the current time (should be set at time of n+1,
  // in the following only (approximations of) rhs's at n+1 are computed
  subSysStatus->updateCurrentTime();

  // reset the achieved flag
  m_data->setAchieved(false);
  m_data->setAchieved(m_data->getNbMaxSteps(0) == 0);

  // update subsystemstatus
  subSysStatus->resetResidual();
  subSysStatus->setFirstStep(true);
  subSysStatus->setLastStep(false);
  
  // We use this convergence status to be able to change the CFL in the pseudo-time
  // iterations in unsteady
  std::auto_ptr<Framework::ConvergenceStatus> cvgst(new Framework::ConvergenceStatus);
  
  // sweeps to solve the nonlinear system
  for(CFuint k = 0; !m_data->isAchieved(); ++k) {
     subSysStatus->setSubIter(k);
     
     *cvgst = subSysStatus->getConvergenceStatus();
    
    // prepare the computation of the residuals
    getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();
    
    // this will make the solvers compute the jacobian only during the first iteration at each time step
    (m_data->freezeJacobian() && k > 1) ? m_data->setDoComputeJacobFlag(false) : m_data->setDoComputeJacobFlag(true);
    
    // this is needed for cases like jacobian free
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag( m_data->getDoComputeJacobFlag() );
    
    // compute the steady space residual
    m_init->execute(); /// @note KVDA: what does this init command do? It's Null by default and it's not present with BDF2...
    
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(0.5);
    if (m_data->getDoComputeJacobFlag()) {
      getConvergenceMethodData()->getCFL()->update(cvgst.get());
    }
    
    // execute intermediate step (add 0.5*rhs^n)
    m_intermediate->execute();
    
    // RHS will be changed here (from steady to pseudoSteady)
    // including the time contribution
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);
    
    // solve the linear system
    getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys),
				  m_data->getNbLSSToSolveAtOnce());
    
    // update the solution
    m_updateSol->execute();

    // synchronize the states and compute the residual norms
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    // Display info over each step of the Newton iterator
    if (m_data->isPrintHistory()) {
      CFout << "Newton Step: " << k+1
           << " L2 dU: " << subSysStatus->getResidual() << "\n";
    }

    // postprocess the solution (typically a limiter)
    getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution();

    // Simple Stop Condition (max number of Steps or L2Norm reached)
    // HERE, if you stop because of residual, then you will do one step too much...
    m_data->setAchieved((k+1 == m_data->getNbMaxSteps(0)) ||
                       (subSysStatus->getResidual() < m_data->getMaxNorm()));

    // set subsystemstatus
    subSysStatus->setLastStep(m_data->isAchieved());
    subSysStatus->setFirstStep(k == 0);
  }

  // back up the spatial rhs
  if(subSysStatus->isSubIterationLastStep()){
    CFLog(VERBOSE, "CrankNicholson::takeStepImpl() => isSubIterationLastStep, SubIter[" 
	   << subSysStatus->getSubIter() << "]\n");
    
    // compute the steady space residual
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(false);
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(0.5);
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(true);
    m_intermediate->execute();
  }
  
  // update mesh
  if(subSysStatus->isMovingMesh()) m_aleUpdate->execute();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

