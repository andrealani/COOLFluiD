#include "Environment/ObjectProvider.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "DiscontGalerkin/DiscontGalerkinSolver.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< DiscontGalerkinSolver,SpaceMethod,DiscontGalerkinModule,1 >
  discontGalerkinSolverProvider("DiscontGalerkinSolver");

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::defineConfigOptions(Config::OptionList& options)
{
	options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
	options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");

	options.addConfigOption< std::vector<std::string> >("BcComds","Types of the boundary conditions commands.");
	options.addConfigOption< std::vector<std::string> >("BcNames","Names for the configuration of the boundary conditions commands.");

	options.addConfigOption< std::string >("SolveFacesCom","Command to solve the problem with Discontinuous Galerkin solver, faces.");
	options.addConfigOption< std::string >("SolveCellsCom","Command to solve the problem with Discontinuous Galerkin solver, cells.");
        options.addConfigOption< std::string >("StabilizationCom","Command to use satbilization in solver.");
	options.addConfigOption< std::string >("SetResidualCom","Command how to set residual.");
	options.addConfigOption< std::string >("SetupCom","Command to initialize Discontinuous Galerkin solver data.");
	options.addConfigOption< std::string >("UnSetupCom","Command to deallocate Discontinuous Galerkin solver data.");
// 	options.addConfigOption< std::string >("SolveTimeDependencies","Command to add time dependet term.");
}

//////////////////////////////////////////////////////////////////////////////

DiscontGalerkinSolver::DiscontGalerkinSolver(const std::string& name) :
  SpaceMethod(name),
  m_setup(),
  m_unsetup(),
  m_solveCells(),
  m_solveFaces(),
  m_setResidual(),
  m_stabilization()
{
  addConfigOptionsTo(this);
  m_data.reset(new DiscontGalerkinSolverData(this));

  cf_assert(m_data.isNotNull());

  // set default value of builder for DiscontGalerkinSolver
  // to be DiscontGalerkinMeshDataBuilder
  m_builder = "DG";
  setParameter( "Builder",   &m_builder );

  // set default global jacobian sparsity
  m_sparsity = "CellCentered";

  m_setupStr   = "StdSetup";
  setParameter( "SetupCom",   &m_setupStr );

  m_unsetupStr = "StdUnSetup";
  setParameter( "UnSetupCom", &m_unsetupStr );

  m_solveCellsStr   = "StdSolveCells";
  setParameter( "SolveCellsCom",   &m_solveCellsStr );

  m_setResidualStr   = "StdComputeResidual";
  setParameter( "SetResidualCom",   &m_setResidualStr );

  m_stabilizationStr   = "Null";
  setParameter( "StabilizationCom",   &m_stabilizationStr );

  m_solveFacesStr   = "Null";
  setParameter( "SolveFacesCom",   &m_solveFacesStr );

  m_initTypeStr = vector<std::string>();
  setParameter("InitComds",&m_initTypeStr);

  m_initNameStr = vector<std::string>();
  setParameter("InitNames",&m_initNameStr);

  m_bcTypeStr = vector<std::string>();
  setParameter("BcComds",&m_bcTypeStr);

  m_bcNameStr = vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);

//   m_timeDependentStr = "StdSolveTimeDependencies";
//   setParameter("SolveTimeDependencies",&m_timeDependentStr);

}

//////////////////////////////////////////////////////////////////////////////

DiscontGalerkinSolver::~DiscontGalerkinSolver()
{
  clearInitComs();
  clearBCComs();
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::setCollaborator( MultiMethodHandle<LinearSystemSolver> lss )
{
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add here configures to the DiscontGalerkinSolver
  configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
    args, m_setup,m_setupStr,m_data );
  configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
    args, m_unsetup,m_unsetupStr,m_data );
  configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
    args, m_solveFaces,m_solveFacesStr,m_data );
  configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
    args, m_solveCells,m_solveCellsStr,m_data );
  configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
    args, m_stabilization,m_stabilizationStr,m_data );
  configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
    args, m_setResidual,m_setResidualStr,m_data );
//   configureCommand< DiscontGalerkinSolverData,DiscontGalerkinSolverCom::PROVIDER >(
//     args, m_timeDependencies,m_timeDependentStr,m_data );

  cf_assert(m_setup.isNotNull());
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_solveCells.isNotNull());
  cf_assert(m_solveFaces.isNotNull());
//   cf_assert(m_timeDependencies.isNotNull());

  clearInitComs();
  clearBCComs();

  cf_assert(m_initTypeStr.size() == m_initNameStr.size());

  m_inits.resize(m_initTypeStr.size());

  for(CFuint i = 0; i < m_inits.size(); ++i) {

    CFLog(INFO, "INIT type = " << m_initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << m_initNameStr[i] << "\n");

    configureCommand<DiscontGalerkinSolverCom,
      DiscontGalerkinSolverData,
      DiscontGalerkinSolverComProvider>
      (args, m_inits[i], m_initTypeStr[i],m_initNameStr[i], m_data);

    cf_assert(m_inits[i].isNotNull());
  }

  cf_assert(m_bcTypeStr.size() == m_bcNameStr.size());

  m_bcs.resize(m_bcTypeStr.size());
  for(CFuint i = 0; i < m_bcs.size(); ++i) {

    CFLog(INFO, "BC type = " << m_bcTypeStr[i] << "\n");
    CFLog(INFO, "BC name = " << m_bcNameStr[i] << "\n");

    configureCommand<DiscontGalerkinSolverCom,
      DiscontGalerkinSolverData,
      DiscontGalerkinSolverComProvider>
      (args, m_bcs[i],m_bcTypeStr[i], m_bcNameStr[i], m_data);

    cf_assert(m_bcs[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::clearInitComs()
{
  vector<SelfRegistPtr<DiscontGalerkinSolverCom> >().swap(m_inits);
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::clearBCComs()
{
  vector<SelfRegistPtr<DiscontGalerkinSolverCom> >().swap(m_bcs);
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::setMethodImpl()
{
  CFAUTOTRACE;

  SpaceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::unsetMethodImpl()
{
  CFAUTOTRACE;

  cf_assert(m_unsetup.isNotNull());
  m_unsetup->execute();
  unsetupCommandsAndStrategies();

  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::initializeSolutionImpl(bool isRestart)
{
   CFAUTOTRACE;

  if (!isRestart) {
     for(CFuint i = 0; i < m_inits.size(); ++i) {
      cf_assert(m_inits[i].isNotNull());
      m_inits[i]->execute();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  Common::Stopwatch<Common::WallTime> stopTimer;
  stopTimer.start();

CFout << "DG begin\n" << CFendl;
if (m_solveFaces.isNotNull()) m_solveFaces->execute();

  cf_assert(m_solveCells.isNotNull());
  m_solveCells->execute();

if (m_stabilization.isNotNull()) m_stabilization->execute();
//   m_timeDependencies->execute();

  applyBC();

  stopTimer.stop();
  CFout << "DG end - DG took " << stopTimer << "s\n" << CFendl;

}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::applyBCImpl()
{
  CFAUTOTRACE;

  for(CFuint i = 0; i < m_bcs.size(); ++i) {
    cf_assert(m_bcs[i].isNotNull());
    m_bcs[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::prepareComputationImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t DiscontGalerkinSolver::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t DiscontGalerkinSolver::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolver::postProcessSolutionImpl()
{
  CFAUTOTRACE;
   cf_assert(m_setResidual.isNotNull());
   m_setResidual->execute();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

