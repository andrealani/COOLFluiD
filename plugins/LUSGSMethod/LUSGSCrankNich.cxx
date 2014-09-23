#include "Common/Stopwatch.hh"

#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUSGSCrankNich.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LUSGSCrankNich,
               ConvergenceMethod,
               LUSGSMethodModule,
               1>
lusgsCrankNichProvider("NonlinearLUSGSCrankNich");

//////////////////////////////////////////////////////////////////////////////

void LUSGSCrankNich::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("BackupPastRhs","Command that backs up the past rhs.");
  options.addConfigOption< std::string >("AddPastRhs","Command that adds the past rhs to the current rhs.");
}

//////////////////////////////////////////////////////////////////////////////

LUSGSCrankNich::LUSGSCrankNich(const std::string& name)
  : LUSGSIterator(name)
{
  addConfigOptionsTo(this);

//   m_setupStr = "CrankNichSetup";
//   setParameter("SetupCom",&m_setupStr);

  m_backupPastRhsStr = "StorePastRhs";
  setParameter("BackupPastRhs",&m_backupPastRhsStr);

  m_addPastRhsStr = "AddPastRhs";
  setParameter("AddPastRhs",&m_addPastRhsStr);
}

//////////////////////////////////////////////////////////////////////////////

LUSGSCrankNich::~LUSGSCrankNich()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSCrankNich::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  LUSGSIterator::configure(args);
  configureNested ( m_data.getPtr(), args );

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>(args, m_backupPastRhs,m_backupPastRhsStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>(args, m_addPastRhs,m_addPastRhsStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSCrankNich::prepare()
{
  CFAUTOTRACE;

  LUSGSIterator::prepare();
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSCrankNich::takeStepImpl()
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

  // STORE THE PAST RHS
  m_data->setForwardSweep(true);
  m_data->setStopSweep(false);
  // Update states set index (it is equal to -1 at this point)
  m_updateStatesSetIndex->execute();
  for (;!m_data->stopSweep();)
  {
    // Compute space residual for the current states set
    m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(0.5);

    // store the past rhs
    m_backupPastRhs->execute();

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

  // SET CURRENT TIME TO TIME AT N+1 (IN THE FOLLOWING ONLY (APPROXIMATIONS OF) THE RHS AT N+1 ARE COMPUTED
  subSysStatus->updateCurrentTime();

  // Compute the block Jacobian matrices at the first time-step and every m_JacobFreezFreq time-steps
  if((subSysStatus->getNbIter() == 1) || ((subSysStatus->getNbIter())%(m_data->getJacobFreezFreq(0))==0))
  {
    timer.restart();

    CFout << "Compute the block Jacobian matrices ...\n";

    // MAKE THE SPACE METHOD COMPUTE THE BLOCK JACOBIAN MATRICES.
    CFLog(VERBOSE,"Start computation of diagonal block Jacobian matrices\n");
    CFLogDebugMed( "LUSGSCrankNich::takeStep(): preparing Computation\n");
    m_data->getCollaborator<SpaceMethod>()->prepareComputation();
    CFLogDebugMed( "LUSGSCrankNich::takeStep(): computing Space Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(0.5);
    CFLogDebugMed( "LUSGSCrankNich::takeStep(): computing the Time Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

    // FACTORIZES THE BLOCK JACOBIAN MATRICES.
    CFLog(VERBOSE,"Starting matrix factorizations\n");
    m_luFactorization->execute();

    // output factorization time
    CFLog(VERBOSE,"Computing and factorizing diagonal block matrices took: " << timer << "s\n");
  }

  CFout << "Starting LU-SGS iterator ...\n";

  // DO THE SGS SWEEPS
  CFuint& k = getConvergenceMethodData()->getConvergenceStatus().iter; // use the convergence status iter of the method
  for(k = 1; !m_data->isAchieved(); ++k)
  {
    timer.restart();

    subSysStatus->setFirstStep( k == 1 );
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());

    // Do forward sweep
    CFLog(VERBOSE,"LUSGSCrankNich::takeStep(): starting forward sweep\n");
    m_data->setForwardSweep(true);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to -1 at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(0.5);

      // add the past rhs
      m_addPastRhs->execute();

      // Compute time residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

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
    CFLog(VERBOSE,"LUSGSCrankNich::takeStep(): starting backward sweep\n");
    m_data->setForwardSweep(false);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to the number of states sets at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(0.5);

      // add the past rhs
      m_addPastRhs->execute();

      // Compute time residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

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
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > LUSGSCrankNich::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result = LUSGSIterator::getStrategyList();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
