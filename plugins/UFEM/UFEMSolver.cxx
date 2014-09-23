#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/UFEMSolver.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< UFEMSolver,SpaceMethod,UFEMPlugin,1 > aUFEMSolverProvider("UFEM");

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
  options.addConfigOption< std::string >("SolveCom","Command to solve the problem with UFEM solver.");
  options.addConfigOption< std::string >("UnSetupCom","Command to deallocate UFEM solver data.");
  options.addConfigOption< std::vector<std::string> >("SetupCom","Command(s) to initialize UFEM solver data.");
  options.addConfigOption< std::string >("CleanCom","Command to clean up before each assembly.");
  options.addConfigOption< std::vector<std::string> >("BcComds","Types of the boundary conditions commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names for the configuration of the boundary conditions commands.");
}

//////////////////////////////////////////////////////////////////////////////

UFEMSolver::UFEMSolver(const std::string& name) :
  SpaceMethod(name),
  m_unsetup(),
  m_solve(),
  m_clean(),
  m_inits(0),
  m_bcs(0)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  m_data.reset(new UFEMSolverData(this));

  cf_assert(m_data.isNotNull());

  // set default value of builder for UFEMSolver
  // to be UFEMMeshDataBuilder
  m_builder = "UFEM";
  // set default global jacobian sparsity
  m_sparsity = "CellVertex";

  m_setupStr.assign(1,"StdSetup");
  setParameter( "SetupCom",   &m_setupStr );

  m_unsetupStr = "StdUnSetup";
  setParameter( "UnSetupCom", &m_unsetupStr );

  m_solveStr   = "StdSolve";
  setParameter( "SolveCom",   &m_solveStr );

  m_cleanStr   = "StdClean";
  setParameter( "CleanCom",   &m_cleanStr );

  m_initTypeStr.clear();
  setParameter("InitComds",&m_initTypeStr);

  m_initNameStr.clear();
  setParameter("InitNames",&m_initNameStr);

  m_bcTypeStr.clear();
  setParameter("BcComds",&m_bcTypeStr);

  m_bcNameStr.clear();
  setParameter("BcNames",&m_bcNameStr);
}

//////////////////////////////////////////////////////////////////////////////

UFEMSolver::~UFEMSolver()
{
  CFAUTOTRACE;
  clearInitComs();
  clearBCComs();
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::clearBCComs()
{
  CFAUTOTRACE;
  m_bcs.clear();
  std::vector<SelfRegistPtr<UFEMSolverCom> >().swap(m_bcs);
}

//////////////////////////////////////////////////////////////////////////////
void UFEMSolver::clearInitComs()
{
  CFAUTOTRACE;
  m_inits.clear();
  std::vector<SelfRegistPtr<UFEMSolverCom> >().swap(m_inits);
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::setCollaborator( MultiMethodHandle<LinearSystemSolver> lss )
{
  CFAUTOTRACE;
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  CFAUTOTRACE;
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // configure setup commands
  m_setup.resize(m_setupStr.size());
  for (CFuint i=0; i<m_setup.size(); ++i) {
    configureCommand< UFEMSolverData,UFEMSolverCom::PROVIDER >(args,m_setup[i],m_setupStr[i],m_data);
    cf_assert(m_setup[i].isNotNull());
  }

  // configure unsetup, solve and clean commands
  configureCommand< UFEMSolverData,UFEMSolverCom::PROVIDER >( args, m_unsetup,m_unsetupStr,m_data );
  configureCommand< UFEMSolverData,UFEMSolverCom::PROVIDER >( args, m_solve,  m_solveStr,  m_data );
  configureCommand< UFEMSolverData,UFEMSolverCom::PROVIDER >( args, m_clean,  m_cleanStr,  m_data );
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_solve.isNotNull());
  cf_assert(m_clean.isNotNull());

  clearInitComs();

  cf_assert(m_initTypeStr.size() == m_initNameStr.size());

  m_inits.resize(m_initTypeStr.size());
  for(CFuint i = 0; i < m_inits.size(); ++i) {

    CFLog(INFO, "INIT type = " << m_initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << m_initNameStr[i] << "\n");

    configureCommand<UFEMSolverCom,
      UFEMSolverData,
      UFEMSolverCom::PROVIDER>( args, m_inits[i],
              m_initTypeStr[i],
              m_initNameStr[i],
              m_data);
    cf_assert(m_inits[i].isNotNull());
  }

  clearBCComs();

  if (m_bcTypeStr.size() != m_bcNameStr.size())
    throw BadValueException(FromHere(), "Number of types of BC does not match number of names for BCs\n");

  m_bcs.resize(m_bcTypeStr.size());
  for(CFuint i = 0; i < m_bcs.size(); ++i) {

    CFLog(INFO, "BC type = " << m_bcTypeStr[i] << "\n");
    CFLog(INFO, "BC name = " << m_bcNameStr[i] << "\n");

    configureCommand<UFEMSolverCom,
      UFEMSolverData,
      UFEMSolverCom::PROVIDER>( args, m_bcs[i],
              m_bcTypeStr[i],
              m_bcNameStr[i],
              m_data);
    cf_assert(m_bcs[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::setMethodImpl()
{
  CFAUTOTRACE;

  SpaceMethod::setMethodImpl();

  setupCommandsAndStrategies();

  for (CFuint i=0; i<m_setup.size(); ++i) {
    cf_assert(m_setup[i].isNotNull());
    m_setup[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::unsetMethodImpl()
{
  CFAUTOTRACE;
  cf_assert(m_unsetup.isNotNull());
  m_unsetup->execute();
  unsetupCommandsAndStrategies();

  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;

  if (!isRestart) // if restart execute all init coms
  {
    for (CFuint i = 0; i < m_inits.size(); ++i)
    {
      cf_assert(m_inits[i].isNotNull());
      m_inits[i]->execute();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
  cf_assert(m_solve.isNotNull());
  m_clean->execute();
  m_solve->execute();
  applyBC();
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::applyBCImpl()
{
  CFAUTOTRACE;
  for(CFuint i = 0; i < m_bcs.size(); ++i)
  {
    cf_assert(m_bcs[i].isNotNull());
    m_bcs[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolver::prepareComputationImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t UFEMSolver::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t UFEMSolver::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace UFEM
}  // namespace COOLFluiD

