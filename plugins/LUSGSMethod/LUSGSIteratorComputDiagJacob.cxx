#include "Common/Stopwatch.hh"

#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUSGSIteratorComputDiagJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LUSGSIteratorComputDiagJacob,
               ConvergenceMethod,
               LUSGSMethodModule,
               1>
lusgsIteratorComputDiagJacobProvider("NonlinearLUSGSIteratorComputDiagJacob");

//////////////////////////////////////////////////////////////////////////////

void LUSGSIteratorComputDiagJacob::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSIteratorComputDiagJacob::LUSGSIteratorComputDiagJacob(const std::string& name)
    : LUSGSIterator(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSIteratorComputDiagJacob::~LUSGSIteratorComputDiagJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIteratorComputDiagJacob::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  LUSGSIterator::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIteratorComputDiagJacob::takeStepImpl()
{
  CFAUTOTRACE;

  // do the prepare step, usually backing up the solution to past states
  prepare();

  // get subsystemstatus
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  // create a timer
  Common::Stopwatch<WallTime> timer;
  Common::Stopwatch<WallTime> total_timer;
  total_timer.restart();

  // Compute the block Jacobian matrices at the first time-step and every m_JacobFreezFreq time-steps
  if((subSysStatus->getNbIter() == 1) || ((subSysStatus->getNbIter())%(m_data->getJacobFreezFreq(0))==0))
  {
    timer.restart();

    CFout << "Compute the block Jacobian matrices ...\n";

    // compute the diagonal block matrices of the Jacobian by perturbation of the states
    CFLog(VERBOSE,"Start computation of diagonal block Jacobian matrices\n");
    m_data->setForwardSweep(true);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to -1 at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(1.0);
      // Compute time residual for the current states set // this contribution should be zero (pastState = state)
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

      // set these booleans to ensure that the correct action is taken in the command that computes the jacobians
      m_data->setStopStatesLoop(true);
      m_data->setStopEqsLoop(true);

      // store the unperturbed residual
      m_diagBlockJacobComputer->execute();

      // loop over the states in the current states set
      m_data->setStopStatesLoop(false);
      for (;!m_data->stopStatesLoop();)
      {
        // loop over the variables in the current state
        m_data->setStopEqsLoop(false);
        for (;!m_data->stopEqsLoop();)
        {
          // call m_diagBlockJacobComputer before computation of perturbed residual
          m_data->setBeforePertResComputation(true);
          m_diagBlockJacobComputer->execute();

          // Compute space residual for the current states set
          m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(1.0);
          // Compute time residual for the current states set
          m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

          // call m_diagBlockJacobComputer before computation of perturbed residual
          m_data->setBeforePertResComputation(false);
          m_diagBlockJacobComputer->execute();
        }
      }

      // Update states set index
      m_updateStatesSetIndex->execute();
    }

    // set the states set index back to -1
    m_data->setForwardSweep(false);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to the number of states sets at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
        // Update states set index
      m_updateStatesSetIndex->execute();
    }

    // Factorizes the block Jacobian matrices.
    CFLog(VERBOSE,"Starting matrix factorizations\n");
    m_luFactorization->execute();

    // output factorization time
    CFLog(VERBOSE,"Computing and factorizing diagonal block matrices took: " << timer << "s\n");
  }

  CFout << "Starting LU-SGS iterator ...\n";

  // Do the SGS sweeps
  CFuint& k = getConvergenceMethodData()->getConvergenceStatus().iter; // use the convergence status iter of the method
  for(k = 1; !m_data->isAchieved(); ++k)
  {

    timer.restart();

    subSysStatus->setFirstStep( k == 1 );
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());

//     m_init->execute();
    // Do forward sweep
    CFLog(VERBOSE,"LUSGSIteratorComputDiagJacob::takeStep(): starting forward sweep\n");
    m_data->setForwardSweep(true);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to -1 at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(1.0);

      // Compute time residual for the current states set
//       m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

      // Compute the solution update for the current states set
      m_computeStatesSetUpdate->execute();

      // Update the solution for the current states set
      m_updateSol->execute();

      // Update states set index
      m_updateStatesSetIndex->execute();
    }

    // Syncronize the states
    ConvergenceMethod::syncGlobalDataComputeResidual(false);

    // reset the residual norms in the local processor
    m_data->getLUSGSNormComputer()->resetResiduals();

    // Do backward sweep
    CFLog(VERBOSE,"LUSGSIteratorComputDiagJacob::takeStep(): starting backward sweep\n");
    m_data->setForwardSweep(false);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to the number of states sets at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(1.0);

      // Compute time residual for the current states set
//       m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

      // Compute the solution update for the current states set
      m_computeStatesSetUpdate->execute();

      // Update the solution for the current states set
      m_updateSol->execute();

      // add contribution of current states set to the residual norms in the local processor
      m_data->getLUSGSNormComputer()->addStatesSetContribution();

      // Update states set index
      m_updateStatesSetIndex->execute();
    }

    // Syncronize the states
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    // output SGS time
    CFLog (VERBOSE, "SGS sweep took: " << timer << "s\n");

    // Gets the residual norm
    getConvergenceMethodData()->getConvergenceStatus().res = subSysStatus->getResidual();

    // Display info over each step of the Nonlinear LU-SGS iterator
    if (m_data->isPrintHistory())
    {
      CFout << "Nonlinear LU-SGS sweep: " << k << " L2 dU: " << getConvergenceMethodData()->getConvergenceStatus().res  << "\n";
    }

//     m_data->setAchieved(m_stopCondControler->isAchieved(getConvergenceMethodData()->getConvergenceStatus()));

    // Simple Stop Condition (max number of Steps or L2Norm reached)
    // HERE, if you stop because of residual, then you will do one step too much...
    m_data->setAchieved((k >= m_data->getNbMaxSweeps(0)) ||
                        (subSysStatus->getResidual() < m_data->getMaxNorm()));

    // set subsystem status
    subSysStatus->setFirstStep(false);
    subSysStatus->setLastStep(m_data->isAchieved());
  }

  CFLog (VERBOSE, "Nonlinear LU-SGS convergence took: " << total_timer << "s\n");

  if (subSysStatus->isSubIterationLastStep()) subSysStatus->updateCurrentTime();
}

std::vector<Common::SafePtr<NumericalStrategy> > LUSGSIteratorComputDiagJacob::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result = LUSGSIterator::getStrategyList();

  return result;
}
//////////////////////////////////////////////////////////////////////////////
    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
