// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "RungeKuttaLS/RungeKuttaLS.hh"
#include "RungeKuttaLS/RKLS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RKLS,
               ConvergenceMethod,
               RungeKuttaLSModule,
               1>
rkConvergenceMethodProvider("RKLS");

//////////////////////////////////////////////////////////////////////////////

void RKLS::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("RungeKuttaStep","Runge-Kutta Step command to run.");
   options.addConfigOption< std::string >("BackupSol","Backup Solution command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

RKLS::RKLS(const std::string& name) :
  ConvergenceMethod(name)
{
   addConfigOptionsTo(this);

   m_data.reset(new RKLSData(this));

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

RKLS::~RKLS()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> RKLS::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> RKLS::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void RKLS::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the RKLSCom's

  configureCommand<RKLSData,RKLSComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<RKLSData,RKLSComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<RKLSData,RKLSComProvider>( args, m_backupSol,m_backupSolStr,m_data);

  configureCommand<RKLSData,RKLSComProvider>( args, m_rungeKuttaStep,m_rungeKuttaStepStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void RKLS::setMethodImpl()
{
  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void RKLS::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void RKLS::takeStepImpl()
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

  // Initialize u0
  m_backupSol->execute();

  // loop over the R-K stages
  for (CFuint i = 0; i < order ; ++i)
  {
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

    // Compute solution of this R-K stage
    m_rungeKuttaStep->execute();


    if (i != order-1)
    {
      // Synchronize the states, do not compute the residual
      ConvergenceMethod::syncGlobalDataComputeResidual(false);

      // postprocess the solution
      m_data->getCollaborator<SpaceMethod>()->postProcessSolution();
    }
  }

  // Synchronize the states, compute the residual
  ConvergenceMethod::syncGlobalDataComputeResidual(true);

  // postprocess the solution
  m_data->getCollaborator<SpaceMethod>()->postProcessSolution();

  // update time to time at end of iteration
  if (dt > 0.)
  {
    SubSystemStatusStack::getActive()->setCurrentTime(m_data->getTimeIterN()+dt);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD
