// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "EmptyIterator.hh"
#include "Environment/ObjectProvider.hh"
#include "EmptyConvergenceMethod/EmptyConvergenceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace EmptyConvergenceMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EmptyIterator,
			    ConvergenceMethod,
			    EmptyConvergenceMethodLib,
			    1>
emptyIteratorMethodProvider("EmptyIterator");

//////////////////////////////////////////////////////////////////////////////

void EmptyIterator::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UpdateSol","Command to update the solution with computed dU.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

EmptyIterator::EmptyIterator(const std::string& name)
  : ConvergenceMethod(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new EmptyIteratorData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "Null";
  setParameter("UnSetupCom",&m_unSetupStr);
  
  m_updateSolStr = "Null";
  setParameter("UpdateSol",&m_updateSolStr);
}

//////////////////////////////////////////////////////////////////////////////

EmptyIterator::~EmptyIterator()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> EmptyIterator::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> EmptyIterator::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void EmptyIterator::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  ConvergenceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add configures to the EmptyIteratorCom's
  
  configureCommand<EmptyIteratorData,EmptyIteratorComProvider>(args, m_setup,m_setupStr,m_data);
  CFLog(VERBOSE, "EmptyIterator::configure() => Command " << m_setupStr << "\n");
  
  configureCommand<EmptyIteratorData,EmptyIteratorComProvider>(args, m_unSetup,m_unSetupStr,m_data);
  CFLog(VERBOSE, "EmptyIterator::configure() => Command " << m_unSetupStr << "\n");
  
  configureCommand<EmptyIteratorData,EmptyIteratorComProvider>(args, m_updateSol,m_updateSolStr,m_data);
  CFLog(VERBOSE, "EmptyIterator::configure() => Command " << m_updateSolStr << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void EmptyIterator::setMethodImpl()
{
  CFAUTOTRACE;

  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();

  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void EmptyIterator::unsetMethodImpl()
{
  m_unSetup->execute();

  unsetupCommandsAndStrategies();

  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void EmptyIterator::takeStepImpl()
{
  CFAUTOTRACE;
  
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  subSysStatus->updateNbIter();
  subSysStatus->updateTimeStep();
  subSysStatus->updateCurrentTime();
  
  CFLog(VERBOSE, "EmptyIterator::takeStepImpl() => iter = " << subSysStatus->getNbIter() << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace EmptyConvergenceMethod

} // namespace COOLFluiD

