// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/Stopwatch.hh"

#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"
#include "Framework/StopConditionController.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/NewtonIterator.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NewtonIterator,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
newtonIteratorConvergenceMethodProvider("NewtonIterator");

//////////////////////////////////////////////////////////////////////////////

void NewtonIterator::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run.");
   options.addConfigOption< std::string >("UpdateSol","Command to update the solution with computed dU.");
   options.addConfigOption< std::string >("InitCom","Command to perform at beggining of each newton step.");
   options.addConfigOption< std::string >("IntermediateCom","Command to perform between the computation of the Space and Time residual.");
   options.addConfigOption< std::string >("PrepareCom","Command to prepare the solution before the iteration process.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run.");
   options.addConfigOption< std::string >("ALEUpdateCom","Command to perform after takeStep when moving mesh.");
}

//////////////////////////////////////////////////////////////////////////////

NewtonIterator::NewtonIterator(const std::string& name)
  : ConvergenceMethod(name)
{
  addConfigOptionsTo(this);

  m_data.reset(new NewtonIteratorData(this));

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

  m_initStr = "Null";
  setParameter("InitCom",&m_initStr);

  m_aleUpdateStr = "Null";
  setParameter("ALEUpdateCom",&m_aleUpdateStr);
 }

//////////////////////////////////////////////////////////////////////////////

NewtonIterator::~NewtonIterator()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NewtonIterator::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> NewtonIterator::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIterator::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ConvergenceMethod::configure(args);
  
  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );
  
  // add configures to the NewtonIteratorCom's

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_prepare,m_prepareStr,m_data);

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_intermediate,m_intermediateStr,m_data);

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_init,m_initStr,m_data);

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_updateSol,m_updateSolStr,m_data);

  configureCommand<NewtonIteratorData,NewtonIteratorComProvider>( args, m_aleUpdate,m_aleUpdateStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIterator::setMethodImpl()
{
  //call the parent
  ConvergenceMethod::setMethodImpl();

  m_data->setLinearSystemSolver(getLinearSystemSolver());
  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIterator::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIterator::prepare()
{
  CFAUTOTRACE;
    
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if (subSysStatus->isSubIterationFirstStep())
  {
    subSysStatus->updateNbIter();
    subSysStatus->updateTimeStep();
//    getConvergenceMethodData()->getCFL()->update();
  }
  
  // prepare to take a time step
  m_prepare->execute();

  getConvergenceMethodData()->getConvergenceStatus().res     = subSysStatus->getResidual();
  getConvergenceMethodData()->getConvergenceStatus().iter    = 0;
  getConvergenceMethodData()->getConvergenceStatus().subiter = 0;
  getConvergenceMethodData()->getConvergenceStatus().time    = subSysStatus->getCurrentTimeDim();

  subSysStatus->setFirstStep(true);
  bool already_achieved = (m_data->getNbMaxSteps(0) == 0) || m_stopCondControler->isAchieved(getConvergenceMethodData()->getConvergenceStatus());
  m_data->setAchieved( already_achieved );
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIterator::takeStepImpl()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "NewtonIterator::takeStepImpl() START\n");

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  
  // AL: temporary solution only checking the solution rate
  const bool solveSystem =
    (subSysStatus->getNbIter() == 0 ||
     (subSysStatus->getNbIter()+1)%getConvergenceMethodData()->getSolvingRate() == 0);
  
  // do the prepare step, usually backing up the solution to past states
  prepare();
  
  Common::Stopwatch<WallTime> timer;
  Common::Stopwatch<WallTime> total_timer;
  total_timer.restart();
  
  // We use this convergence status to be able to change the CFL in the pseudo-time
  // iterations in unsteady
  std::auto_ptr<Framework::ConvergenceStatus> cvgst(new Framework::ConvergenceStatus);
  // use the convergence status iter of the method
  CFuint& k = getConvergenceMethodData()->getConvergenceStatus().iter;
  
  if (getConvergenceMethodData()->onlyPreprocessSolution()) {
    getMethodData()->getCollaborator<SpaceMethod>()->setOnlyPreprocessSolution(true);
    // start by preprocessing the solution (once for the whole simulation) w/o solving anything
    getMethodData()->getCollaborator<SpaceMethod>()->preProcessSolution();
  }
  
  if (solveSystem) {
    CFLog(VERBOSE, "NewtonIterator::takeStepImpl() for " << getName() << "\n");
    
  for(k = 1; !m_data->isAchieved(); ++k)
  {    
    *cvgst = subSysStatus->getConvergenceStatus();
	
    timer.restart();
    
    subSysStatus->setFirstStep( k == 1 );
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());
    
    m_init->execute();
    CFLog(VERBOSE, "NewtonIterator::takeStep(): preparing Computation\n");
    getMethodData()->getCollaborator<SpaceMethod>()->prepareComputation();
   
    CFLog(VERBOSE, "NewtonIterator::takeStep(): m_data->freezeJacobian() " << m_data->freezeJacobian() << "\n");
    // this will make the solvers compute the jacobian only during the first iteration at each time step
    (m_data->freezeJacobian() && k > 1) ? m_data->setDoComputeJacobFlag(false) : m_data->setDoComputeJacobFlag(true);
    
    // this is needed for cases like jacobian free
    getMethodData()->getCollaborator<SpaceMethod>()->setComputeJacobianFlag( m_data->getDoComputeJacobFlag() );
    
    CFLog(VERBOSE, "NewtonIterator::takeStep(): computing Space Residual\n");
    // compute the steady space residual
    getMethodData()->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
    
    CFLog(VERBOSE, "NewtonIterator::takeStep(): before second update CFL\n");
    
    if (m_data->getDoComputeJacobFlag()) {
      getConvergenceMethodData()->getCFL()->update(cvgst.get());
    }
    
    // do an intermediate step, useful for some special types of temporal discretization
    CFLog(VERBOSE, "NewtonIterator::takeStep(): calling Intermediate step\n");
    m_intermediate->execute();
    
    CFLog(VERBOSE, "NewtonIterator::takeStep(): computing the Time Residual\n");
    getMethodData()->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

    CFLog(VERBOSE, "Assembling linear system took: " << timer << "s\n");

    // RHS will be changed here (from steady to pseudoSteady)
    // including the time contribution
    CFLog(VERBOSE, "NewtonIterator::takeStep(): solving the linear system\n");
    timer.restart();

    // solve the linear system
    getLinearSystemSolver().apply(mem_fun(&LinearSystemSolver::solveSys), 
				  m_data->getNbLSSToSolveAtOnce());
    
    CFLog(VERBOSE, "Solving linear system took: " << timer << "s\n");

    /// @todo each processor should print in separate files
    if(m_data->isSaveSystemToFile())
    {
      std::string prefix = "system-";
      std::string suffix =  "." + StringOps::to_str(k) + ".dat";
      for (CFuint i = 0; i < getLinearSystemSolver().size(); ++i)
      {
        getLinearSystemSolver()[i]->printToFile(prefix,suffix);
      }
    }

    CFLog(VERBOSE, "NewtonIterator::takeStep(): updating the solution\n");
    m_updateSol->execute();
    
    // synchronize the states and compute the residual norms
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    getMethodData()->getCollaborator<SpaceMethod>()->postProcessSolution();
    getConvergenceMethodData()->getConvergenceStatus().res = subSysStatus->getResidual();

    // Display info over each step of the Newton iterator
    if (m_data->isPrintHistory())
    {
      CFout << "Newton Step: " << k << " L2 dU: " << getConvergenceMethodData()->getConvergenceStatus().res  << " CFL: " << getConvergenceMethodData()->getCFL()->getCFLValue() << "\n";
    }

    m_data->setAchieved(m_stopCondControler->isAchieved(getConvergenceMethodData()->getConvergenceStatus()));
  }
  
  CFLog (VERBOSE, "Newton convergence took: " << total_timer << "s\n");
  }
  
  if (subSysStatus->isSubIterationLastStep()) subSysStatus->updateCurrentTime();
  
  CFLog(VERBOSE, "NewtonIterator::takeStepImpl() END\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

