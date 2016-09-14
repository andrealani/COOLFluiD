// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "BDF2.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "Framework/CFL.hh"
#include "Framework/SpaceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<BDF2,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
bdf2MethodProvider("BDF2");

//////////////////////////////////////////////////////////////////////////////

BDF2::BDF2(const std::string& name)
  : NewtonIterator(name)
{
  // new default commands are set for the BDF2 scheme
  m_setupStr   = "BDF2Setup";
  m_unSetupStr = "CrankNichUnSetup";
  m_intermediateStr = "BDF2Intermediate";
}

//////////////////////////////////////////////////////////////////////////////

BDF2::~BDF2()
{
}

//////////////////////////////////////////////////////////////////////////////

void BDF2::prepare()
{
  CFLog(VERBOSE, "BDF2::prepare() START\n");
  
  // update subsystemstatus if first subiteration step
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if (subSysStatus->isSubIterationFirstStep())
  {
    subSysStatus->updateNbIter();
    subSysStatus->updateTimeStep();
    // getConvergenceMethodData()->getCFL()->update();
  }
  
  // prepare to take a time step
  m_prepare->execute();

  // reset the achieved flag
  m_data->setAchieved(false);
  m_data->setAchieved(m_data->getNbMaxSteps(0) == 0);

  // update subsystemstatus
  subSysStatus->resetResidual();
  subSysStatus->setFirstStep(true);
  subSysStatus->setLastStep(false);
  
  CFLog(VERBOSE, "BDF2::prepare() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void BDF2::takeStepImpl()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "BDF2::takeStepImpl() START\n");
  
  // do the prepare step, usually backing up the solution to past states
  prepare ();

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  // get time steps
  const CFreal newDT = subSysStatus->getDT();
  const CFreal oldDT = subSysStatus->getPreviousDT();

  // compute factors related to BDF2 scheme
  const CFreal alpha = oldDT / newDT;
  CFreal xi = 1. / ( 1.+ alpha);
  CFreal theta = 1.;

  // for the first step, backward Euler is used
  if(subSysStatus->getNbIter() == 1)
  {
    xi = 0.;
    theta = 1.0;
  }

  std::auto_ptr<Framework::ConvergenceStatus> cvgst(new Framework::ConvergenceStatus);
  // sweeps to solve the nonlinear system
  for(CFuint k = 0; !m_data->isAchieved(); ++k)
  {
    *cvgst = subSysStatus->getConvergenceStatus();
    
    CFLog(VERBOSE, "BDF2::takeStep(): preparing Computation\n");
    
    // prepare the computation of the residuals
    getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();
   
    CFLog(VERBOSE, "BDF2::takeStep(): m_data->freezeJacobian() " << m_data->freezeJacobian() << "\n");
    // this will make the solvers compute the jacobian only during the first iteration at each time step
    (m_data->freezeJacobian() && k > 1) ? m_data->setDoComputeJacobFlag(false) : m_data->setDoComputeJacobFlag(true);
    
    // this is needed for cases like jacobian free or in case the jacobian needs to be frozen for k>=1
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(m_data->getDoComputeJacobFlag());
    
    CFLog(VERBOSE, "BDF2::takeStep(): m_data->getDoComputeJacobFlag() " << m_data->getDoComputeJacobFlag() << ", k = " << k << "\n");
    CFLog(VERBOSE, "BDF2::takeStep(): computing Space Residual\n");
    // compute the steady space residual
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(theta);
    
    if ( m_data->getDoComputeJacobFlag()) {
      getConvergenceMethodData()->getCFL()->update(cvgst.get());
    }
    
    CFLog(VERBOSE, "BDF2::takeStep(): computing space residual norm\n");
    
    // The norm of the space residual will only be computed if the flag 
    // "OutputSpaceResidual" is set to true in ConvergenceMethod.cxx
    getConvergenceMethodData()->computeSpaceResidualNorm();

    CFLog(VERBOSE, "BDF2::takeStep(): computing Time Residual\n");

    // RHS will be changed here (from steady to pseudoSteady)
    // including the time contribution
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(xi);
    CFLog(VERBOSE, "BDF2::takeStep(): solving the linear system\n");
    
    // solve the linear system
    getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys),
				  m_data->getNbLSSToSolveAtOnce());
    
    CFLog(VERBOSE, "BDF2::takeStep(): updating the solution\n");

    CFLog(VERBOSE, "BDF2::takeStep(): updating the solution\n");
    m_updateSol->execute();
    
    // synchronize the states and compute the residual norms
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    // Display info over each step of the Newton iterator
    if (m_data->isPrintHistory())
    {
      CFout << "BDF2 Step: " << k+1 << " L2 dU: " << subSysStatus->getResidual() << "\n";
    }
    
    CFLog(VERBOSE, "BDF2::takeStep(): post process the solution\n"); 
    
    // postprocess the solution (typically a limiter)
    getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution();

    // Simple Stop Condition (max number of Steps or L2Norm reached)
    // HERE, if you stop because of residual, then you will do one step too much...
    m_data->setAchieved((k+1 == m_data->getNbMaxSteps(0)) ||
			(subSysStatus->getResidual() < m_data->getMaxNorm()));

    // set subsystemstatus
    subSysStatus->setFirstStep(false);
    subSysStatus->setLastStep(m_data->isAchieved());
  }

  // back up the time contribution
  if(subSysStatus->isSubIterationLastStep())
  {
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(false);
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(xi);
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag(true);
    subSysStatus->updateCurrentTime();
  }

  // update mesh
  if(subSysStatus->isMovingMesh()) m_aleUpdate->execute();

  CFLog(VERBOSE, "BDF2::takeStepImpl() END\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

