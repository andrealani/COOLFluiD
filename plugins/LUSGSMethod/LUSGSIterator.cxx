#include "Common/Stopwatch.hh"

#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUSGSIterator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LUSGSIterator,
               ConvergenceMethod,
               LUSGSMethodModule,
               1>
lusgsIteratorProvider("NonlinearLUSGSIterator");

//////////////////////////////////////////////////////////////////////////////

void LUSGSIterator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SetupCom","SetupCommand to run.");
  options.addConfigOption< std::string >("UpdateSol","Command to update the solution with computed dU.");
  options.addConfigOption< std::string >("InitCom","Command to perform at beggining of each (nonlinear) LU-SGS step.");
  options.addConfigOption< std::string >("IntermediateCom","Command to perform between the computation of the Space and Time residual.");
  options.addConfigOption< std::string >("PrepareCom","Command to prepare the solution before the iteration process.");
  options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run.");
  options.addConfigOption< std::string >("ALEUpdateCom","Command to perform after takeStep when moving mesh.");
  options.addConfigOption< std::string >("UpdateStatesSetIdx","Command that updates the states set index.");
  options.addConfigOption< std::string >("LUFactorization","Command to perform the LU factorization.");
  options.addConfigOption< std::string >("ComputeSolUpdate","Command to solve the two triangular systems after the LU factorization.");
  options.addConfigOption< std::string >("ComputeJacobians","Command for the computation of the diagonal block Jacobians.");
}

//////////////////////////////////////////////////////////////////////////////

LUSGSIterator::LUSGSIterator(const std::string& name)
  : ConvergenceMethod(name)
{
   addConfigOptionsTo(this);

  m_data.reset(new LUSGSIteratorData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_prepareStr = "StdPrepare";
  setParameter("PrepareCom",&m_prepareStr);

  m_updateSolStr = "UpdateStatesSetSolution";
  setParameter("UpdateSol",&m_updateSolStr);

  m_intermediateStr = "Null";
  setParameter("IntermediateCom",&m_intermediateStr);

  m_initStr = "Null";
  setParameter("InitCom",&m_initStr);

  m_aleUpdateStr = "Null";
  setParameter("ALEUpdateCom",&m_aleUpdateStr);

  m_updateStatesSetIndexStr = "UpdateStatesSetIndex";
  setParameter("UpdateStatesSetIdx",&m_updateStatesSetIndexStr);

  m_luFactorizationStr = "LUFactPivot";
  setParameter("LUFactorization",&m_luFactorizationStr);

  m_computeStatesSetUpdateStr = "ComputeStatesSetUpdatePivot";
  setParameter("ComputeSolUpdate",&m_computeStatesSetUpdateStr);

  m_diagBlockJacobComputerStr = "DiagBlockJacobMatrByPert";
  setParameter("ComputeJacobians",&m_diagBlockJacobComputerStr);
}

//////////////////////////////////////////////////////////////////////////////

LUSGSIterator::~LUSGSIterator()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> LUSGSIterator::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> LUSGSIterator::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIterator::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ConvergenceMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add configures to the LUSGSIteratorCom's

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_prepare,m_prepareStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_intermediate,m_intermediateStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_init,m_initStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_updateSol,m_updateSolStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_aleUpdate,m_aleUpdateStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_updateStatesSetIndex,m_updateStatesSetIndexStr,m_data);

  CFout << "Using " << m_luFactorizationStr << "\n";
  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_luFactorization,m_luFactorizationStr,m_data);

  CFout << "Using " << m_computeStatesSetUpdateStr << "\n";
  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_computeStatesSetUpdate,m_computeStatesSetUpdateStr,m_data);

  configureCommand<LUSGSIteratorData,LUSGSIteratorComProvider>( args, m_diagBlockJacobComputer,m_diagBlockJacobComputerStr,m_data);

}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIterator::setMethodImpl()
{
  //call the parent
  ConvergenceMethod::setMethodImpl();

//   m_data->setLinearSystemSolver(getLinearSystemSolver());
  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIterator::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();
  //call the parent
  ConvergenceMethod::unsetMethodImpl();
}


//////////////////////////////////////////////////////////////////////////////

void LUSGSIterator::prepare()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if (subSysStatus->isSubIterationFirstStep())
  {
    subSysStatus->updateNbIter();
    subSysStatus->updateTimeStep();
    getConvergenceMethodData()->getCFL()->update();
  }

  m_prepare->execute();

  getConvergenceMethodData()->getConvergenceStatus().res     = subSysStatus->getResidual();
  getConvergenceMethodData()->getConvergenceStatus().iter    = 0;
  getConvergenceMethodData()->getConvergenceStatus().subiter = 0;
  getConvergenceMethodData()->getConvergenceStatus().time    = subSysStatus->getCurrentTimeDim();

  SubSystemStatusStack::getActive()->setFirstStep(true);
  SubSystemStatusStack::getActive()->setLastStep(false);
  bool already_achieved = (m_data->getNbMaxSweeps(0) == 0) || m_stopCondControler->isAchieved(getConvergenceMethodData()->getConvergenceStatus());
  m_data->setAchieved( already_achieved );
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSIterator::takeStepImpl()
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

    // Makes the space method compute the block Jacobian matrices.
    CFLog(VERBOSE,"Start computation of diagonal block Jacobian matrices\n");
    CFLogDebugMed( "LUSGSIterator::takeStep(): preparing Computation\n");
    m_data->getCollaborator<SpaceMethod>()->prepareComputation();
    CFLogDebugMed( "LUSGSIterator::takeStep(): computing Space Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
    //  fill in socket_diagBlockJacobMatr

    CFLogDebugMed( "LUSGSIterator::takeStep(): computing the Time Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);
    //  update socket_diagBlockJacobMatr

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
    CFLog(VERBOSE,"LUSGSIterator::takeStep(): starting forward sweep\n");
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
    CFLog(VERBOSE,"LUSGSIterator::takeStep(): starting backward sweep\n");
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

    // Display info over each sweep of the Nonlinear LU-SGS iterator
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

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > LUSGSIterator::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result;

  // add strategies here
  result.push_back(m_data->getNormComputer().d_castTo<NumericalStrategy>());

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
