// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RK2.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SpaceMethod.hh"
#include "RungeKutta2/RungeKutta2.hh"
#include "Framework/CFL.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RK2,
         ConvergenceMethod,
         RungeKutta2Module,
         1>
rk2ConvergenceMethodProvider("RK2");

//////////////////////////////////////////////////////////////////////////////

void RK2::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("PredictorStep","Predictor Step command to run.");
   options.addConfigOption< std::string >("BackupSol","Backup Solution command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("CorrectorStep","Corrector Step command to run.");
}

//////////////////////////////////////////////////////////////////////////////

RK2::RK2(const std::string& name) :
  ConvergenceMethod(name)
{
   addConfigOptionsTo(this);

   m_data.reset(new RK2Data(this));

   m_setupStr = "StdSetup";
   setParameter("SetupCom",&m_setupStr);



   m_unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&m_unSetupStr);



   m_predictorStepStr = "PredictorStep";
   setParameter("PredictorStep",&m_predictorStepStr);



   m_backupSolStr = "StdBackupSol";
   setParameter("BackupSol",&m_backupSolStr);



   m_correctorStepStr = "CorrectorStep";
   setParameter("CorrectorStep",&m_correctorStepStr);


}

//////////////////////////////////////////////////////////////////////////////

RK2::~RK2()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> RK2::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> RK2::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void RK2::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the RK2Com's

  configureCommand<RK2Data,RK2ComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<RK2Data,RK2ComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<RK2Data,RK2ComProvider>( args, m_backupSol,m_backupSolStr,m_data);

  configureCommand<RK2Data,RK2ComProvider>( args, m_predictorStep,m_predictorStepStr,m_data);

  configureCommand<RK2Data,RK2ComProvider>( args, m_correctorStep,m_correctorStepStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void RK2::setMethodImpl()
{
  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void RK2::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void RK2::takeStepImpl()
{
  CFAUTOTRACE;

  SubSystemStatusStack::getActive()->updateNbIter();
  SubSystemStatusStack::getActive()->updateTimeStep();
  getConvergenceMethodData()->getCFL()->update();

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  m_backupSol->execute();

  m_data->getCollaborator<SpaceMethod>()->prepareComputation();
  m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
  m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

  m_predictorStep->execute();

  ConvergenceMethod::syncGlobalDataComputeResidual(false);

  m_data->getCollaborator<SpaceMethod>()->postProcessSolution();

  if (dt > 0.)
  {
    SubSystemStatusStack::getActive()->
        setCurrentTime(SubSystemStatusStack::getActive()->getCurrentTime() + 0.5*dt);
  }

  m_data->getCollaborator<SpaceMethod>()->prepareComputation();
  m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
  m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

  m_correctorStep->execute();

  ConvergenceMethod::syncGlobalDataComputeResidual(true);

  m_data->getCollaborator<SpaceMethod>()->postProcessSolution();

  if (dt > 0.)
  {
    SubSystemStatusStack::getActive()->
        setCurrentTime(SubSystemStatusStack::getActive()->getCurrentTime() + 0.5*dt);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD
