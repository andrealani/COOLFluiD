#include "Framework/VarRegistry.hh"

#include "Environment/ObjectProvider.hh"

#include "Framework/SpaceMethod.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/CFL.hh"

#include "RKRD/RKRD.hh"
#include "RKRD/RungeKuttaRD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RKRD,
               ConvergenceMethod,
               RKRDModule,
               1>
RKRD_Provider("RKRD");

//////////////////////////////////////////////////////////////////////////////

void RKRD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("BackupCom","Backup Command to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("UpdateCom","Update Command to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("ShiftCom","Shift Command to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

RKRD::RKRD(const std::string& name) :
  ConvergenceMethod(name)
{
   addConfigOptionsTo(this);

   m_data.reset(new RKRDData(this));

   m_setupStr = "StdSetup";
   setParameter("SetupCom",&m_setupStr);

   m_backupStr = "Backup";
   setParameter("BackupCom",&m_backupStr);

   m_updateStr = "Update";
   setParameter("UpdateCom",&m_updateStr);

   m_shiftStr = "Shift";
   setParameter("ShiftCom",&m_shiftStr);

   m_unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&m_unSetupStr);
}

//////////////////////////////////////////////////////////////////////////////

RKRD::~RKRD()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> RKRD::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> RKRD::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void RKRD::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the RKRDCom's

  configureCommand<RKRDData,RKRDComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<RKRDData,RKRDComProvider>( args, m_backup,m_backupStr,m_data);

  configureCommand<RKRDData,RKRDComProvider>( args, m_update,m_updateStr,m_data);

  configureCommand<RKRDData,RKRDComProvider>( args, m_shift,m_shiftStr,m_data);

  configureCommand<RKRDData,RKRDComProvider>( args, m_unSetup,m_unSetupStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void RKRD::setMethodImpl()
{
  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  m_setup->execute();

  Common::SafePtr<SubSystemStatus> subsys_status = SubSystemStatusStack::getActive();
  CFuint * kstep = new CFuint();
  subsys_status->getVarRegistry()->registVar<CFuint>("kstep", kstep );
}

//////////////////////////////////////////////////////////////////////////////

void RKRD::unsetMethodImpl()
{
  Common::SafePtr<SubSystemStatus> subsys_status = SubSystemStatusStack::getActive();
  CFuint * kstep = subsys_status->getVarRegistry()->unregistVar<CFuint>("kstep");
  deletePtr(kstep);

  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void RKRD::takeStepImpl()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> subsys_status = SubSystemStatusStack::getActive();

  // update number of iterations, time step and CFL
  subsys_status->updateNbIter();
  subsys_status->updateTimeStep();
  getConvergenceMethodData()->getCFL()->update();

  // get time step at this iteration
  const CFreal cur_time = SubSystemStatusStack::getActive()->getCurrentTime();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  // Get the order of the method
  const CFuint order = m_data->getOrder() ;

  // copy the states U_n to all U_k
  m_backup->execute();

  // loop over the R-K stages
  CFuint& k = m_data->K();
  for ( k = 0; k < order ; ++k )
  {
    // set current RK step
    subsys_status->getVarRegistry()->setVar<CFuint>("kstep",k);

    // Compute the RHS
    Common::SafePtr<SpaceMethod> rd = m_data->getCollaborator<SpaceMethod>();
    rd->prepareComputation();
    rd->computeSpaceResidual(1.0);
    rd->computeTimeResidual(1.0);

    // intermidiate update and shift
    m_shift->execute();

    // synchronize parallel data and compute residual
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    // post process solution
    rd->postProcessSolution();
  }

  // final update of the states
  m_update->execute();

  // update time to time at end of iteration
  if (dt > 0.)
  {
    SubSystemStatusStack::getActive()->setCurrentTime(cur_time + dt);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD
