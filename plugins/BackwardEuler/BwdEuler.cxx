// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/Stopwatch.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/SpaceMethod.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"

#include "BackwardEuler/BackwardEuler.hh"
#include "BackwardEuler/BwdEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<BwdEuler,
               ConvergenceMethod,
               BackwardEulerModule,
               1>
bwdEulerConvergenceMethodProvider("BwdEuler");

//////////////////////////////////////////////////////////////////////////////

void BwdEuler::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UpdateSol","Update Solution command.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

BwdEuler::BwdEuler(const std::string& name) :
  ConvergenceMethod(name),
  m_data(new BwdEulerData(this)),
  m_setup(),
  m_unSetup(),
  m_updateSol()
{
  addConfigOptionsTo(this);
  m_setupStr = "StdSetup";
   setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_updateSolStr = "UpdateSol";
  setParameter("UpdateSol",&m_updateSolStr);
}

//////////////////////////////////////////////////////////////////////////////

BwdEuler::~BwdEuler()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> BwdEuler::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> BwdEuler::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void BwdEuler::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add configures to the BwdEulerCom's
  configureCommand<BwdEulerData,BwdEulerComProvider>(args, m_setup, m_setupStr,m_data);

  configureCommand<BwdEulerData,BwdEulerComProvider>(args, m_unSetup, m_unSetupStr,m_data);

  configureCommand<BwdEulerData,BwdEulerComProvider>(args, m_updateSol, m_updateSolStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void BwdEuler::setMethodImpl()
{
  //call the parent
  ConvergenceMethod::setMethodImpl();

  m_data->setLinearSystemSolver(getLinearSystemSolver());
  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void BwdEuler::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void BwdEuler::takeStepImpl()
{
  CFAUTOTRACE;

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep())
  {
    SubSystemStatusStack::getActive()->updateNbIter();
    SubSystemStatusStack::getActive()->updateTimeStep();
    // getConvergenceMethodData()->getCFL()->update();
    
    getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();
  }
  
  Common::Stopwatch<Common::WallTime> timer;
  timer.restart();

  getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
  getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

  timer.stop();
  CFLog(INFO, "BwdEuler : assembled linear system in " << timer << "s\n");
  timer.restart();

  getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys),
				m_data->getNbLSSToSolveAtOnce());
  timer.stop();
  CFLog(INFO, "BwdEuler : solved linear system in " << timer << "s\n");

  m_updateSol->execute();

  ConvergenceMethod::syncGlobalDataComputeResidual(true);

  getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution();

  if(SubSystemStatusStack::getActive()->isSubIterationLastStep())
  {
    SubSystemStatusStack::getActive()->updateCurrentTime();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

