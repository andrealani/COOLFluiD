// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "RungeKutta/RungeKutta.hh"
#include "RungeKutta/RK.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RK,
               ConvergenceMethod,
               RungeKuttaModule,
               1>
rkConvergenceMethodProvider("RK");

//////////////////////////////////////////////////////////////////////////////

void RK::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("RungeKuttaStep","RungeKutta Step command to run.");
   options.addConfigOption< std::string >("BackupSol","Backup Solution command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

RK::RK(const std::string& name) :
  ConvergenceMethod(name)
{
   addConfigOptionsTo(this);

   m_data.reset(new RKData(this));

   m_setupStr = "StdSetup";
   setParameter("SetupCom",&m_setupStr);

   m_unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&m_unSetupStr);

   m_rungeKuttaStepStr = "RungeKuttaStep";
   setParameter("RungeKuttaStep",&m_rungeKuttaStepStr);

   m_backupSolStr = "StdBackupSol";
   setParameter("BackupSol",&m_backupSolStr);
}

//////////////////////////////////////////////////////////////////////////////

RK::~RK()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> RK::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> RK::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void RK::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the RKCom's

  configureCommand<RKData,RKComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<RKData,RKComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<RKData,RKComProvider>( args, m_backupSol,m_backupSolStr,m_data);

  configureCommand<RKData,RKComProvider>( args, m_rungeKuttaStep,m_rungeKuttaStepStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void RK::setMethodImpl()
{
  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void RK::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void RK::takeStepImpl()
{
  CFAUTOTRACE;

  // update number of iterations, time step and CFL
  SubSystemStatusStack::getActive()->updateNbIter();
  SubSystemStatusStack::getActive()->updateTimeStep();
  getConvergenceMethodData()->getCFL()->update();

  // store time at beginning of iteration
  m_data->setTimeIterN(SubSystemStatusStack::getActive()->getCurrentTime());

  // get time step at this iteration
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  // Get the order of the method
  const CFuint order = m_data->getOrder() ;

  // Initialize U0 and tempu
  m_data->setIsLastStep(false);
  m_data->setIsFirstStep(true);
  m_backupSol->execute();
  m_data->setIsFirstStep(false);

  // loop over the R-K stages
  for (CFuint i = 0; i < order ; ++i) {

    // set current step in the method data
    m_data->setCurrentStep(i);

    // update the time for the current stage
    if (dt > 0.)
    {
      SubSystemStatusStack::getActive()->setCurrentTime(m_data->getTimeIterN()+m_data->getGamma(i)*dt);
    }

    // Compute the RHS
    m_data->getCollaborator<SpaceMethod>()->prepareComputation();
    m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
    m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

    // Compute Ui and Un+1(i)
    m_rungeKuttaStep->execute();
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    m_data->getCollaborator<SpaceMethod>()->postProcessSolution();
  }

  // Update the states
  m_data->setIsLastStep(true);
  m_backupSol->execute();

  // update time to time at end of iteration
  if (dt > 0.)
  {
    SubSystemStatusStack::getActive()->setCurrentTime(m_data->getTimeIterN()+dt);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD
