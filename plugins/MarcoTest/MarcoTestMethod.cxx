// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MarcoTest/MarcoTestMethod.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SubSystemStatus.hh"
#include "MarcoTest/MarcoTest.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MarcoTestMethod,
               ConvergenceMethod,
               MarcoTestModule,
               1>
fwdEulerConvergenceMethodProvider("MarcoTestMethod");

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethod::defineConfigOptions(Config::OptionList& options)
{
  // Options that you can set up in the input file 
  options.addConfigOption< std::string >("SetupCom"," Command that creates data structure");
  options.addConfigOption< std::string >("AlgoCom"," Command that implement shocking");
  options.addConfigOption< std::string >("UnSetupCom"," Command that unsetup the data");
}

//////////////////////////////////////////////////////////////////////////////

MarcoTestMethod::MarcoTestMethod(const std::string& name)
  : ConvergenceMethod(name),
    m_data(), m_setup(), m_algo(), m_unsetup()
{
  addConfigOptionsTo(this);
  m_data.reset(new MarcoTestMethodData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);
  
  m_algoStr = "Shock";
  setParameter("AlgoCom",&m_algoStr); 
  
  m_unsetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unsetupStr);
}
      
//////////////////////////////////////////////////////////////////////////////

MarcoTestMethod::~MarcoTestMethod()
{
  CFLog(VERBOSE, "MarcoTestMethod::~MarcoTestMethod()\n");
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> MarcoTestMethod::getMethodData () const
{
   return m_data.getPtr();
 }

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> MarcoTestMethod::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  ConvergenceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add configures to the MarcoTestMethodCom's
  
  configureCommand<MarcoTestMethodData,MarcoTestMethodComProvider>
    (args, m_setup,m_setupStr,m_data);
  CFLog(VERBOSE, "MarcoTestMethod::configure() => Command " << m_setupStr << "\n");
  
  configureCommand<MarcoTestMethodData,MarcoTestMethodComProvider>
    (args, m_algo,m_algoStr,m_data);
  CFLog(VERBOSE, "MarcoTestMethod::configure() => Command " << m_algoStr << "\n");

  configureCommand<MarcoTestMethodData,MarcoTestMethodComProvider>
    (args, m_unsetup,m_unsetupStr,m_data);
  CFLog(VERBOSE, "MarcoTestMethod::configure() => Command " << m_unsetupStr << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethod::setMethodImpl()
{
  CFAUTOTRACE;

  ConvergenceMethod::setMethodImpl();
  
  m_setup->execute(); 
 
  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethod::unsetMethodImpl()
{ 
  //  unsetupCommandsAndStrategies();
  unsetupCommandsAndStrategies();

  m_unsetup->execute();

  ConvergenceMethod::unsetMethodImpl(); 
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethod::takeStepImpl()
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  
  subSysStatus->updateNbIter();
  subSysStatus->updateTimeStep();
  
  // // Commands have execute, setup and unsetup Marco
  m_algo->execute();
  ConvergenceMethod::syncGlobalDataComputeResidual(false);
  subSysStatus->updateCurrentTime();
  
  // m_data->getNodalStatesExtrapolator()->extrapolateInAllNodes();
}
    
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > MarcoTestMethod::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result;
  result.push_back(m_data->getGeoDataComputer().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getNodalStatesExtrapolator().d_castTo<NumericalStrategy>());
  return result;
}
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest

} // namespace COOLFluiD

