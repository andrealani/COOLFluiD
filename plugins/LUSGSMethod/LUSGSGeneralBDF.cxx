#include "Common/Stopwatch.hh"
#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUSGSGeneralBDF.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LUSGSGeneralBDF,
               ConvergenceMethod,
               LUSGSMethodModule,
               1>
lusgsGeneralBDFProvider("NonlinearLUSGSGeneralBDF");

//////////////////////////////////////////////////////////////////////////////

LUSGSGeneralBDF::LUSGSGeneralBDF(const std::string& name)
  : LUSGSIterator(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSGeneralBDF::~LUSGSGeneralBDF()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSGeneralBDF::prepare()
{
  CFAUTOTRACE;

  LUSGSIterator::prepare();
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSGeneralBDF::takeStepImpl()
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

  // set coefficients to 1. The coefficients used for the general BDF scheme with variable time step are computed in the space method
  const CFreal xi = 1.0;
  const CFreal theta = 1.0;

  // SET CURRENT TIME TO TIME AT N+1 (IN THE FOLLOWING ONLY (APPROXIMATIONS OF) THE RHS AT N+1 ARE COMPUTED
  subSysStatus->updateCurrentTime();

// Compute the block Jacobian matrices at the first time-step and every m_JacobFreezFreq time-steps
  if((subSysStatus->getNbIter() == 1) || ((subSysStatus->getNbIter())%(m_data->getJacobFreezFreq(0))==0))
  {
    timer.restart();

    CFout << "Compute the block Jacobian matrices ...\n";

    // MAKE THE SPACE METHOD COMPUTE THE BLOCK JACOBIAN MATRICES.
    CFLog(VERBOSE,"Start computation of diagonal block Jacobian matrices\n");
    CFLogDebugMed( "LUSGSGeneralBDF::takeStep(): preparing Computation\n");
    m_data->getCollaborator<SpaceMethod>()->prepareComputation();
    CFLogDebugMed( "LUSGSGeneralBDF::takeStep(): computing Space Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(theta);
    CFLogDebugMed( "LUSGSGeneralBDF::takeStep(): computing the Time Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(xi);

    // FACTORIZES THE BLOCK JACOBIAN MATRICES.
    CFLog(VERBOSE,"Starting matrix factorizations\n");
    m_luFactorization->execute();

    // output factorization time
    CFLog(VERBOSE,"Computing and factorizing diagonal block matrices took: " << timer << "s\n");

  }

  CFout << "Starting LU-SGS iterator ...\n";

  // DO THE SGS SWEEPS
  CFuint k = 0;
  for(; k < m_data->getNbMaxSweeps(0) && !m_data->isAchieved(); ++k)
  {
    timer.restart();

    // set subsystem status
    subSysStatus->setFirstStep( k == 0 );
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());

    // Do forward sweep
    CFLog(VERBOSE,"LUSGSGeneralBDF::takeStep(): starting forward sweep\n");
    m_data->setForwardSweep(true);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to -1 at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(theta);

      // Compute time residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(xi);

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
    CFLog(VERBOSE,"LUSGSGeneralBDF::takeStep(): starting backward sweep\n");
    m_data->setForwardSweep(false);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to the number of states sets at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(theta);

      // Compute time residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(xi);

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
      CFout << "Nonlinear LU-SGS sweep: " << k+1 << " L2 dU: " << getConvergenceMethodData()->getConvergenceStatus().res  << "\n";
    }

    // Simple Stop Condition (L2Norm reached)
    // HERE, if you stop because of residual, then you will do one step too much...
    m_data->setAchieved(subSysStatus->getResidual() < m_data->getMaxNorm());

    // set subsystem status
    subSysStatus->setFirstStep(false);
    subSysStatus->setLastStep(m_data->isAchieved());
  }

  // BACK UP THE TIME CONTRIBUTION
  m_data->setForwardSweep(true);
  m_data->setStopSweep(false);
  // Update states set index (it is equal to -1 at this point)
  m_updateStatesSetIndex->execute();
  for (;!m_data->stopSweep();)
  {
    // back up the time rhs for the current states set
    m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

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

  CFLog (VERBOSE, "Nonlinear LU-SGS convergence took: " << total_timer << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
