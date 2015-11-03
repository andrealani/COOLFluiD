// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "FwdEuler.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SpaceMethod.hh"
#include "Common/CFLog.hh"
#include "Framework/SubSystemStatus.hh"
#include "ForwardEuler/ForwardEuler.hh"
#include "Framework/CFL.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FwdEuler,
               ConvergenceMethod,
               ForwardEulerLib,
               1>
fwdEulerConvergenceMethodProvider("FwdEuler");

//////////////////////////////////////////////////////////////////////////////

void FwdEuler::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UpdateSol","Command to update the solution with computed dU.");
   options.addConfigOption< std::string >("IntermediateCom","Command to perform between the computation of the Space and Time residual.");
   options.addConfigOption< std::string >("PrepareCom","Command to prepare the solution before the iteration process.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

FwdEuler::FwdEuler(const std::string& name)
  : ConvergenceMethod(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new FwdEulerData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_prepareStr = "StdPrepare";
  setParameter("PrepareCom",&m_prepareStr);

  m_updateSolStr = "StdUpdateSol";
  setParameter("UpdateSol",&m_updateSolStr);

  m_intermediateStr = "Null";
  setParameter("IntermediateCom",&m_intermediateStr);
}

//////////////////////////////////////////////////////////////////////////////

FwdEuler::~FwdEuler()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> FwdEuler::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> FwdEuler::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void FwdEuler::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  ConvergenceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add configures to the FwdEulerCom's

  configureCommand<FwdEulerData,FwdEulerComProvider>(args, m_setup,m_setupStr,m_data);
  CFLog(VERBOSE, "FwdEuler::configure() => Command " << m_setupStr << "\n");
  
  configureCommand<FwdEulerData,FwdEulerComProvider>(args, m_unSetup,m_unSetupStr,m_data);
  CFLog(VERBOSE, "FwdEuler::configure() => Command " << m_unSetupStr << "\n");

  configureCommand<FwdEulerData,FwdEulerComProvider>(args, m_prepare,m_prepareStr,m_data);
  CFLog(VERBOSE, "FwdEuler::configure() => Command " << m_prepareStr << "\n");

  configureCommand<FwdEulerData,FwdEulerComProvider>(args, m_intermediate,m_intermediateStr,m_data);
  CFLog(VERBOSE, "FwdEuler::configure() => Command " << m_intermediateStr << "\n");

  configureCommand<FwdEulerData,FwdEulerComProvider>(args, m_updateSol,m_updateSolStr,m_data);
  CFLog(VERBOSE, "FwdEuler::configure() => Command " << m_updateSolStr << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void FwdEuler::setMethodImpl()
{
  CFAUTOTRACE;

  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();

  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FwdEuler::unsetMethodImpl()
{
  m_unSetup->execute();

  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void FwdEuler::takeStepImpl()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  subSysStatus->updateNbIter();
  subSysStatus->updateTimeStep();
  getConvergenceMethodData()->getCFL()->update();

  subSysStatus->setFirstStep(true);

  // do a prepare step, usually backing up the solution to pastStates
  CFLog(VERBOSE, "ForwardEuler::takeStep(): calling Prepare step\n");
  if (m_prepare->isNotNull()) { m_prepare->execute(); }

  getConvergenceMethodData()->getConvergenceStatus().res     = subSysStatus->getResidual();
  getConvergenceMethodData()->getConvergenceStatus().iter    = 0;
  getConvergenceMethodData()->getConvergenceStatus().subiter = 0;
  getConvergenceMethodData()->getConvergenceStatus().time    = subSysStatus->getCurrentTimeDim();

  bool already_achieved = m_stopCondControler->isAchieved(getConvergenceMethodData()->getConvergenceStatus());
  m_data->setAchieved( already_achieved );

  CFuint& k = getConvergenceMethodData()->getConvergenceStatus().iter; // use the convergence status iter of the method
  for(k = 1; !m_data->isAchieved(); ++k)
  {
    if (k>1) { subSysStatus->setFirstStep(false); }

    // raise the maximum DT, the space method will compute the maximum allowed and lower it
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());

    CFLog(VERBOSE, "ForwardEuler::takeStep(): computing Space Residual\n");
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual (1.0);

    // do an intermediate step, useful for some special
    // types of temporal discretization
    CFLog(VERBOSE, "ForwardEuler::takeStep(): calling Intermediate step\n");
    m_intermediate->execute();

    CFLog(VERBOSE, "ForwardEuler::takeStep(): computing the Time Residual\n");
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

    CFLog(VERBOSE, "ForwardEuler::takeStep(): updating the solution\n");

    if (m_data->getDoUpdateSolution()) {
      m_updateSol->execute();
    }

    CFLog(VERBOSE, "ForwardEuler::syncGlobalDataComputeResidual()\n");

    ConvergenceMethod::syncGlobalDataComputeResidual(true);
 
    CFLog(VERBOSE, "getCollaborator<SpaceMethod>()->postProcessSolution()\n");

    getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution();

    getConvergenceMethodData()->getConvergenceStatus().res = subSysStatus->getResidual();

    //RealVector l2_norm = m_data->getNorm();
    // getConvergenceMethodData()->getConvergenceStatus().res = l2_norm.emax();

    // Display info over each step of the Forward Euler
    if (m_data->isPrintHistory())
    {
      CFout << "ForwardEuler Step: " << k << " L2 dU: " << getConvergenceMethodData()->getConvergenceStatus().res << "\n";
    }

    m_data->setAchieved(m_stopCondControler->isAchieved(getConvergenceMethodData()->getConvergenceStatus()));
  }

  ConvergenceMethod::syncGlobalDataComputeResidual(false);
  subSysStatus->updateCurrentTime();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler

} // namespace COOLFluiD

