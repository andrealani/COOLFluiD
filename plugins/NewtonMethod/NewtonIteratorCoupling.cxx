// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/Stopwatch.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/CFL.hh"
#include "Framework/SpaceMethod.hh"

#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/NewtonIteratorCoupling.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NewtonIteratorCoupling,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
newtonIteratorCouplingMethodProvider("NewtonIteratorCoupling");

//////////////////////////////////////////////////////////////////////////////

NewtonIteratorCoupling::NewtonIteratorCoupling(const std::string& name)
  : NewtonIterator(name)
{
  m_updateSolStr = "UpdateSolCoupling";
}

//////////////////////////////////////////////////////////////////////////////

NewtonIteratorCoupling::~NewtonIteratorCoupling()
{
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIteratorCoupling::takeStepImpl()
{
  CFAUTOTRACE;

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    SubSystemStatusStack::getActive()->updateNbIter();
    SubSystemStatusStack::getActive()->updateTimeStep();
    getConvergenceMethodData()->getCFL()->update();

  }
  // do a prepare step, usually backing up the solution to
  // pastStates
  CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): calling Prepare step\n");
  m_prepare->execute();

  m_data->setAchieved(false);
  SubSystemStatusStack::getActive()->setFirstStep(true);

  MultiMethodHandle<LinearSystemSolver> lss = m_data->getLinearSystemSolver();
  const CFuint nbLSS = lss.size();

  // if all the linear systems need 0 iterations then consider
  // achieved the convergence
  CFuint count = 0;
  for (CFuint iLSS = 0; iLSS < nbLSS; ++iLSS) {
    if(m_data->getNbMaxSteps(iLSS) == 0) {
      count++;
    }
    else {
      break;
    }
  }
  if (count == nbLSS) {
    m_data->setAchieved(true);
  }

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  Common::Stopwatch<WallTime> timer;

  if (!m_data->isAchieved()) {

    // solve and update one system at a time
    for (CFuint iLSS = 0; iLSS < nbLSS; ++iLSS) {
      // set the data describing the current equation subsystem
      const vector<CFuint>& currEqs = *lss[iLSS]->getEquationIDs();
      const CFuint nbEqsLSS = currEqs.size();
      const CFuint start = currEqs[0];

      // set the equation subsystem descriptor
      PhysicalModelStack::getActive()->
        setEquationSubSysDescriptor(start, nbEqsLSS, iLSS);

      // reset to not achieved before considering the new system
      m_data->setAchieved(false);

      for(CFuint k = 0; !m_data->isAchieved(); ++k) {
    	timer.restart();

      if (k>0) { subSysStatus->setFirstStep(false); }

      subSysStatus->setMaxDT(10.0e10);

      m_init->execute();

      CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): preparing Computation\n");
      getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();

      CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): computing Space Residual\n");
      getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);

      // do an intermediate step, useful for some special
      // types of temporal discretization
      CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): calling Intermediate step\n");
      m_intermediate->execute();

      CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): computing the Time Residual\n");
      getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

      timer.stop();
      CFLogDebugMin("Jacobian matrix computation took: " << timer << "s\n");

      CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): solving the linear system\n");
      timer.restart();

      lss[iLSS]->solveSys();

      timer.stop();
      CFLogDebugMin("Jacobian matrix solution took: " << timer << "s\n");

	    /// @TODO each processor should print in separate files
	    if(m_data->isSaveSystemToFile())
	    {
        std::string prefix = "system-";
        std::string suffix =  "-" + StringOps::to_str(k) + ".dat";
        lss[iLSS]->printToFile(prefix,suffix);
	    }

      CFLogDebugMed( "NewtonIteratorCoupling::takeStep(): updating the solution\n");
      m_updateSol->execute();

      /// @TODO the synchronization may be correctly placed here !!!!
	    ConvergenceMethod::syncGlobalDataComputeResidual(true);

//       getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution));

	    // Display info over each step of the Newton iterator
	    if (m_data->isPrintHistory()) {
	      CFout << "Newton Step: " << k+1
		    << " L2 dU: " << subSysStatus->getResidual() << "\n";
	    }

	    // Simple Stop Condition (max number of Steps or L2Norm reached
	    m_data->setAchieved((k+1 == m_data->getNbMaxSteps(iLSS)));
	    /// @TODO check for the residual is not possible for the moment because
	      // we should check different variables residuals, one for each
	      // subsystem ...
	      //  ||  (subSysStatus->getResidual() < m_data->getMaxNorm()));
       }

    } // end loop iLSS
  }

  if(SubSystemStatusStack::getActive()->isSubIterationLastStep()){
    SubSystemStatusStack::getActive()->updateCurrentTime();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD


