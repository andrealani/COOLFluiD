// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "BDF2_InitCN.hh"
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

Environment::ObjectProvider<BDF2_InitCN,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
bdf2_InitCNProvider("BDF2_InitCN");

//////////////////////////////////////////////////////////////////////////////

BDF2_InitCN::BDF2_InitCN(const std::string& name) :
    NewtonIterator(name)
{
  // new default commands are set for the BDF2scheme
  m_setupStr   = "BDF2Setup";
  m_unSetupStr = "CrankNichUnSetup";
  m_intermediateStr = "BDF2Intermediate";

  m_prepare1stStepStr      = "CN1stStepPrepare";
  m_intermediate1stStepStr = "CN1stStepIntermediate";
}

//////////////////////////////////////////////////////////////////////////////

BDF2_InitCN::~BDF2_InitCN()
{
}

//////////////////////////////////////////////////////////////////////////////

void BDF2_InitCN::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  NewtonIterator::configure(args);

  // add configures to the NewtonIteratorCom's
  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_prepare1stStep     ,m_prepare1stStepStr     ,m_data);
  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_intermediate1stStep,m_intermediate1stStepStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void BDF2_InitCN::takeStepImpl()
{
  // update subsystemstatus if first subiteration step
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  if(subSysStatus->isSubIterationFirstStep())
  {
    subSysStatus->updateNbIter();
    subSysStatus->updateTimeStep();
    getConvergenceMethodData()->getCFL()->update();
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

  //for the first step, CranckNicholson is used
  if(subSysStatus->getNbIter() == 1)
  {
    // prepare for the CranckNicholson step (backup rhs^0)
    getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(0.5);
    m_prepare1stStep->execute();

    // sweeps to solve the nonlinear system
    for(CFuint k = 0; !m_data->isAchieved(); ++k)
    {
      // prepare the computation of the residuals
      getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();

      // compute the steady space residual
      getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(0.5);

      // add rhs^0
      m_intermediate1stStep->execute();

      // RHS will be changed here (from steady to pseudoSteady)
      // including the time contribution
      getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(0.0);

      // solve the linear system
      getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys));

      // update the solution
      m_updateSol->execute();

      // synchronize the states and compute the residual norms
      ConvergenceMethod::syncGlobalDataComputeResidual(true);

      // Display info over each step of the Newton iterator
      if (m_data->isPrintHistory())
      {
        CFout << "BDF2 Step: 1 (Cranck-Nicholson is used) L2 dU: "
              << subSysStatus->getResidual() << "\n";
      }

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
  }
  else // regular BDF2_InitCN
  {
    // get time steps
    const CFreal newDT = subSysStatus->getDT();
    const CFreal oldDT = subSysStatus->getPreviousDT();

    // compute factors related to BDF2_InitCN scheme
    const CFreal alpha = oldDT / newDT;
    CFreal xi = 1. / ( 1.+ alpha);
    CFreal theta = 1.;

    // sweeps to solve the nonlinear system
    for(CFuint k = 0; !m_data->isAchieved(); ++k)
    {
      // prepare the computation of the residuals
      getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();

      // compute the steady space residual
      getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(theta);

      // RHS will be changed here (from steady to pseudoSteady)
      // including the time contribution
      getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(xi);

      // solve the linear system
      getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys));

      // update the solution
      m_updateSol->execute();

      // synchronize the states and compute the residual norms
      ConvergenceMethod::syncGlobalDataComputeResidual(true);

      // Display info over each step of the Newton iterator
      if (m_data->isPrintHistory())
      {
        CFout << "BDF2 Step: " << k+1
              << " L2 dU: " << subSysStatus->getResidual() << "\n";
      }

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
  }

  if(subSysStatus->isSubIterationLastStep())
  {
    // back up the time contribution
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(0.0);
    subSysStatus->updateCurrentTime();
  }

  if(subSysStatus->isMovingMesh()) m_aleUpdate->execute();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

