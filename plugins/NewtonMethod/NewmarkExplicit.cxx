// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewmarkExplicit.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "Framework/CFL.hh"
#include "Framework/SpaceMethod.hh"
#include "Common/CFLog.hh"
#include "Common/Stopwatch.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NewmarkExplicit,
                            ConvergenceMethod,
                            NewtonMethodModule,
                            1>
NewmarkExplicitMethodProvider("NewmarkExplicit");

//////////////////////////////////////////////////////////////////////////////

NewmarkExplicit::NewmarkExplicit(const std::string& name)
  : NewtonIterator(name)
{
  // new default commands are set for the NewmarkExplicit scheme
  m_initStr = "ResetSystem";
  m_setupStr   = "NewmarkSetup";
  m_unSetupStr = "NewmarkUnSetup";
  m_prepareStr = "NewmarkPrepare";
  m_updateSolStr = "NewmarkExplicitUpdateSol";
}

//////////////////////////////////////////////////////////////////////////////

NewmarkExplicit::~NewmarkExplicit()
{
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkExplicit::takeStepImpl()
{
  CFAUTOTRACE;

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    SubSystemStatusStack::getActive()->updateNbIter();
    SubSystemStatusStack::getActive()->updateTimeStep();
    getConvergenceMethodData()->getCFL()->update();
  }

  // do a prepare step, usually backing up the solution to
  // pastStates
  CFLogDebugMed( "NewmarkExplicit::takeStep(): calling Prepare step\n");
  m_prepare->execute();

  m_data->setAchieved(false);
  SubSystemStatusStack::getActive()->setFirstStep(true);
  if(m_data->getNbMaxSteps(0) == 0) m_data->setAchieved(true);
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  Common::Stopwatch<WallTime> timer;

  for(CFuint k = 0; !m_data->isAchieved(); ++k) {

    timer.restart();

    if (k>0) subSysStatus->setFirstStep(false);
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());

    m_init->execute();
    CFLogDebugMed( "NewtonIterator::takeStep(): preparing Computation\n");

	getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();

    CFLogDebugMed( "NewtonIterator::takeStep(): computing Space Residual\n");
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);

    // do an intermediate step, useful for some special
    // types of temporal discretization
    CFLogDebugMed( "NewtonIterator::takeStep(): calling Intermediate step\n");
    m_intermediate->execute();

    CFLogDebugMed( "NewtonIterator::takeStep(): computing the Time Residual\n");
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);


    timer.stop();
    CFLogDebugMin("Jacobian matrix computation took: " << timer << "s\n");

    CFLogDebugMed( "NewtonIterator::takeStep(): solving the linear system\n");
    timer.restart();
    getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys));
    timer.stop();
    CFLogDebugMin("Jacobian matrix solution took: " << timer << "s\n");

    /// @TODO each processor should print in separate files
    if(m_data->isSaveSystemToFile()) {
      std::string prefix = "system-";
      std::string suffix =  "." + StringOps::to_str(k) + ".dat";
      for (CFuint i = 0; i < getLinearSystemSolver().size(); ++i) {
             getLinearSystemSolver()[i]->printToFile(prefix,suffix);
      }
    }

    CFLogDebugMed( "NewtonIterator::takeStep(): updating the solution\n");
    m_updateSol->execute();

    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution();

    // Display info over each step of the Newton iterator
    if (m_data->isPrintHistory()) {
      CFout << "Newton Step: " << k+1
           << " L2 dU: " << subSysStatus->getResidual() << "\n";
    }

    // Simple Stop Condition (max number of Steps or L2Norm reached
    m_data->setAchieved((k+1 == m_data->getNbMaxSteps(0)) ||
                       (subSysStatus->getResidual() < m_data->getMaxNorm()));
    /// @TODO it is dangerous to compare just with one residual
      // this is not fully correct if two systems are solved sequentially
      // residuals of more than one variables may be checked instead
  }

  if(SubSystemStatusStack::getActive()->isSubIterationLastStep()){
    SubSystemStatusStack::getActive()->updateCurrentTime();
  }
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

