// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SAMGLSS/SAMGLSS.hh"
#include "SAMGLSS/SAMGLSSModule.hh"
#include "Framework/BlockAccumulator.hh"
#include "Environment/ObjectProvider.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace SAMGLSS {

Environment::ObjectProvider< SAMGLSS,LinearSystemSolver,SAMGLSSModule,1 >
  samglssLSSMethodProvider("SAMGLSS");

//////////////////////////////////////////////////////////////////////////////

void SAMGLSS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >( "SetupCom",
    "Setup Command to run. This command seldomly needs overriding." );
  options.addConfigOption< std::string >( "UnSetupCom",
    "UnSetup Command to run. This command seldomly needs overriding." );
  options.addConfigOption< std::string >( "SysSolver",
    "Command that solves the linear system." );
}

//////////////////////////////////////////////////////////////////////////////


SAMGLSS::SAMGLSS(const std::string& name) :
  LinearSystemSolver(name)
{
  CFAUTOTRACE;

  m_data.reset(new SAMGLSSData(getMaskArray(),getNbSysEquations(),this ));
  cf_assert(m_data.getPtr() != CFNULL);
  
  addConfigOptionsTo(this);

  m_setupStr = "StdSetup";
  m_solveSysStr = "StdSolveSys";
  m_unSetupStr = "StdUnSetup";
  setParameter("SetupCom",&m_setupStr);
  setParameter("SysSolver",&m_solveSysStr);
  setParameter("UnSetupCom",&m_unSetupStr);
}

//////////////////////////////////////////////////////////////////////////////

SAMGLSS::~SAMGLSS()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSS::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  LinearSystemSolver::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the SAMGLSS
  configureCommand< SAMGLSSData,SAMGLSSComProvider >(
    m_setup,m_setupStr,m_data );
  configureCommand< SAMGLSSData,SAMGLSSComProvider >(
    m_unSetup,m_unSetupStr,m_data );
  configureCommand< SAMGLSSData,SAMGLSSComProvider >(
    m_solveSys,m_solveSysStr,m_data );
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSS::solveSysImpl()
{
  CFAUTOTRACE;
  cf_assert(isSetup());
  cf_assert(isConfigured());
  m_solveSys->execute();
}

//////////////////////////////////////////////////////////////////////////////

BlockAccumulator* SAMGLSS::createBlockAccumulator(
  const CFuint nbRows, const CFuint nbCols, const CFuint subBlockSize ) const
{
  CFAUTOTRACE;
  return new BlockAccumulator(
    nbRows, nbCols, subBlockSize, m_lssData->getLocalToGlobalMapping() );
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSS::printToFile(const std::string prefix, const std::string suffix)
{
  CFAUTOTRACE;
  cf_assert(isSetup());
  cf_assert(isConfigured());
  m_data.getPtr()->printToFile(prefix,suffix);
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSS::setMethodImpl()
{
  CFAUTOTRACE;

  LinearSystemSolver::setMethodImpl();

  m_setup->setup();
  m_setup->execute();

  m_solveSys->setup();
  m_unSetup->setup();

}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSS::unsetMethodImpl()
{
  CFAUTOTRACE;
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  LinearSystemSolver::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr< Framework::MethodData > SAMGLSS::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

