// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Pardiso/Pardiso.hh"
#include "Pardiso/PardisoModule.hh"
#include "Framework/BlockAccumulator.hh"
#include "Environment/ObjectProvider.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Pardiso {

Environment::ObjectProvider< Pardiso,LinearSystemSolver,PardisoModule,1 >
  pardisoLSSMethodProvider("Pardiso");

//////////////////////////////////////////////////////////////////////////////

void Pardiso::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >( "SetupCom",   "Setup Command to run. This command seldomly needs overriding." );
  options.addConfigOption< std::string >( "UnSetupCom", "UnSetup Command to run. This command seldomly needs overriding." );
  options.addConfigOption< std::string >( "SysSolver",  "Command that solves the linear system." );
}

//////////////////////////////////////////////////////////////////////////////


Pardiso::Pardiso(const std::string& name) :
  LinearSystemSolver(name)
{
  CFAUTOTRACE;

  m_data.reset(new PardisoData(getMaskArray(),getNbSysEquations(),this ));
  cf_assert(m_data.getPtr() != CFNULL);
  
  addConfigOptionsTo(this);

  m_setupStr    = "StdSetup";
  m_solveSysStr = "StdSolveSys";
  m_unSetupStr  = "StdUnSetup";
  setParameter("SetupCom",&m_setupStr);
  setParameter("SysSolver",&m_solveSysStr);
  setParameter("UnSetupCom",&m_unSetupStr);
}

//////////////////////////////////////////////////////////////////////////////

Pardiso::~Pardiso()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Pardiso::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  LinearSystemSolver::configure(args);
  configureNested(m_data.getPtr(),args);

  // add here configures
  configureCommand< PardisoData,PardisoComProvider >(args,m_setup,m_setupStr,m_data);
  configureCommand< PardisoData,PardisoComProvider >(args,m_unSetup,m_unSetupStr,m_data);
  configureCommand< PardisoData,PardisoComProvider >(args,m_solveSys,m_solveSysStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void Pardiso::solveSysImpl()
{
  CFAUTOTRACE;
  cf_assert(isSetup());
  cf_assert(isConfigured());
  m_solveSys->execute();
}

//////////////////////////////////////////////////////////////////////////////

BlockAccumulator* Pardiso::createBlockAccumulator(
  const CFuint nbRows, const CFuint nbCols, const CFuint subBlockSize ) const
{
  CFAUTOTRACE;
  return new BlockAccumulator(
    nbRows, nbCols, subBlockSize, m_lssData->getLocalToGlobalMapping() );
}

//////////////////////////////////////////////////////////////////////////////

void Pardiso::printToFile(const std::string prefix, const std::string suffix)
{
  CFAUTOTRACE;
  cf_assert(isSetup());
  cf_assert(isConfigured());
  m_data.getPtr()->printToFile(prefix,suffix);
}

//////////////////////////////////////////////////////////////////////////////

void Pardiso::setMethodImpl()
{
  CFAUTOTRACE;

  LinearSystemSolver::setMethodImpl();

  m_setup->setup();
  m_setup->execute();

  m_solveSys->setup();
  m_unSetup->setup();

}

//////////////////////////////////////////////////////////////////////////////

void Pardiso::unsetMethodImpl()
{
  CFAUTOTRACE;
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  LinearSystemSolver::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr< Framework::MethodData > Pardiso::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Pardiso
} // namespace COOLFluiD

