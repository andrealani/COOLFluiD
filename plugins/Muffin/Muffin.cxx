
#include "Environment/ObjectProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Muffin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< Muffin,SpaceMethod,MuffinModule,1 > camusProvider("Muffin");

//////////////////////////////////////////////////////////////////////////////

Muffin::Muffin(const std::string& name) : SpaceMethod(name)
{
  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  // create data
  m_data.reset(new MuffinData(this));
  cf_assert(m_data.isNotNull());

  // set mesh data builder and global jacobian sparsity
  m_builder  = "Muffin";
  m_sparsity = "CellVertex";

  // (standard commands)
  m_linkStr    = "StdLink";
  m_setupStr   = "StdSetup";
  m_solveStr   = "StdSolve";
  m_unsetupStr = "StdUnSetup";
  setParameter("LinkCom",&m_linkStr);
  setParameter("SetupCom",&m_setupStr);
  setParameter("SolveCom",&m_solveStr);
  setParameter("UnSetupCom",&m_unsetupStr);

  // (special commands)
  m_specialcomds.assign(1,"ComSave");
  setParameter("SpecialComds",&m_specialcomds);

  // (loops, systems, coupling and boundary conditions commands)
  m_loopcomds.clear();
  m_loopnames.clear();
  m_syscomds.clear();
  m_sysnames.clear();
  m_cccomds.clear();
  m_ccnames.clear();
  m_bccomds.clear();
  m_bcnames.clear();
  setParameter("LoopComds",&m_loopcomds);
  setParameter("LoopNames",&m_loopnames);
  setParameter("SystemComds",&m_syscomds);
  setParameter("SystemNames",&m_sysnames);
  setParameter("CcComds",&m_cccomds);
  setParameter("CcNames",&m_ccnames);
  setParameter("BcComds",&m_bccomds);
  setParameter("BcNames",&m_bcnames);
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::defineConfigOptions(Config::OptionList& options)
{
  // (standard commands)
  options.addConfigOption< std::string >("LinkCom","Standard commmand to link Systems to coupling and boundary conditions (default \"StdLink\")");
  options.addConfigOption< std::string >("SetupCom","Standard commmand to build method required data (default \"StdSetup\")");
  options.addConfigOption< std::string >("SolveCom","Standard commmand to solve testcase (default \"StdSolve\")");
  options.addConfigOption< std::string >("UnSetupCom","Standard commmand to dispose of built data (default \"StdUnSetup\")");

  // (special commands)
  options.addConfigOption< std::vector< std::string > >("SpecialComds","Special commands that can be called within loops (default < \"ComSave\" >)");

  // (loops, systems, coupling and boundary conditions commands)
  options.addConfigOption< std::vector< std::string > >("LoopComds","Types of loop commands");
  options.addConfigOption< std::vector< std::string > >("LoopNames","Names of loop commands");
  options.addConfigOption< std::vector< std::string > >("SystemComds","Types of systems commands");
  options.addConfigOption< std::vector< std::string > >("SystemNames","Names of systems commands");
  options.addConfigOption< std::vector< std::string > >("CcComds","Types of coupling conditions commands");
  options.addConfigOption< std::vector< std::string > >("CcNames","Names of coupling conditions commands");
  options.addConfigOption< std::vector< std::string > >("BcComds","Types of boundary conditions commands");
  options.addConfigOption< std::vector< std::string > >("BcNames","Names of boundary conditions commands");
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  SpaceMethod::configure(args);

  configureNested(m_data.getPtr(),args);

  // configure standard commands
  configureCommand< MuffinData,MuffinCom::PROVIDER >(args,m_link,m_linkStr,m_data);
  configureCommand< MuffinData,MuffinCom::PROVIDER >(args,m_setup,m_setupStr,m_data);
  configureCommand< MuffinData,MuffinCom::PROVIDER >(args,m_solve,m_solveStr,m_data);
  configureCommand< MuffinData,MuffinCom::PROVIDER >(args,m_unsetup,m_unsetupStr,m_data);
  cf_assert(m_link.isNotNull());
  cf_assert(m_setup.isNotNull());
  cf_assert(m_solve.isNotNull());
  cf_assert(m_unsetup.isNotNull());

  // configure special commands
  m_vcomm_specials.resize(m_specialcomds.size());
  for (CFuint i=0; i<m_specialcomds.size(); ++i) {
    m_data->log("Special = " + m_specialcomds[i]);
    configureCommand< MuffinData,MuffinCom::PROVIDER >(
      args,m_vcomm_specials[i],m_specialcomds[i],m_data );
    cf_assert(m_vcomm_specials[i].isNotNull());
  }

  // configure loop commands
  if (m_loopcomds.size()!=m_loopnames.size())
    m_data->err("LoopComds and LoopNames should have the same size");
  m_vcomm_loops.resize(m_loopcomds.size());
  for (CFuint i=0; i<m_loopcomds.size(); ++i) {
    m_data->log("Loop = " + m_loopnames[i] + " (" + m_loopcomds[i] + ")");
    configureCommand< MuffinData,MuffinCom::PROVIDER >(
      args,m_vcomm_loops[i],m_loopcomds[i],m_loopnames[i],m_data );
    cf_assert(m_vcomm_loops[i].isNotNull());
  }

  // configure systems commands
  if (m_syscomds.size()!=m_sysnames.size())
    m_data->err("SystemComds and SystemNames should have the same size");
  m_vcomm_sys.resize(m_syscomds.size());
  for (CFuint i=0; i<m_syscomds.size(); ++i) {
    m_data->log("System = " + m_sysnames[i] + " (" + m_syscomds[i] + ")");
    configureCommand< MuffinData,MuffinCom::PROVIDER >(
      args,m_vcomm_sys[i],m_syscomds[i],m_sysnames[i],m_data );
    cf_assert(m_vcomm_sys[i].isNotNull());
  }

  // configure coupling conditions commands
  if (m_cccomds.size()!=m_ccnames.size())
    m_data->err("CcComds and CcNames should have the same size");
  m_vcomm_cc.resize(m_cccomds.size());
  for (CFuint i=0; i<m_cccomds.size(); ++i) {
    m_data->log("CC = " + m_ccnames[i] + " (" + m_cccomds[i] + ")");
    configureCommand< MuffinData,MuffinCom::PROVIDER >(
      args,m_vcomm_cc[i],m_cccomds[i],m_ccnames[i],m_data );
    cf_assert(m_vcomm_cc[i].isNotNull());
  }

  // configure boundary conditions commands
  if (m_bccomds.size()!=m_bcnames.size())
    m_data->err("BcComds and BcNames should have the same size");
  m_vcomm_bc.resize(m_bccomds.size());
  for (CFuint i=0; i<m_bccomds.size(); ++i) {
    m_data->log("BC = " + m_bcnames[i] + " (" + m_bccomds[i] + ")");
    configureCommand< MuffinData,MuffinCom::PROVIDER >(
      args,m_vcomm_bc[i],m_bccomds[i],m_bcnames[i],m_data );
    cf_assert(m_vcomm_bc[i].isNotNull());
  }

  // set special, loop, system, coupling and boundary condition commands in data
  m_data->setCommands(m_vcomm_specials,m_vcomm_loops,m_vcomm_sys,m_vcomm_cc,m_vcomm_bc);
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::setMethodImpl()
{
  SpaceMethod::setMethodImpl();

  m_link->execute();
  m_setup->execute();
  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::unsetMethodImpl()
{
  CFAUTOTRACE;
  m_unsetup->execute();
  unsetupCommandsAndStrategies();

  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
  m_solve->execute();
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::applyBCImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Muffin::prepareComputationImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t Muffin::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t Muffin::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

