#include "FiniteElementMethod.hh"
#include "Environment/ObjectProvider.hh"
#include "ComputeResidualStrategy.hh"
#include "ComputeJacobStrategy.hh"
#include "ConvectiveEntity.hh"
#include "DiffusiveEntity.hh"
#include "InertiaEntity.hh"
#include "ComputeConvectiveTerm.hh"
#include "ComputeDiffusiveTerm.hh"
#include "ComputeLinearSourceTerm.hh"
#include "ComputeIndepSourceTerm.hh"
#include "ComputeInertiaTerm.hh"

#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FiniteElementMethod,
         SpaceMethod,
         FiniteElementModule,
         1>
finiteElementMethodProvider("FiniteElementMethod");

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
   options.addConfigOption< std::string >("ComputeSpaceResidual","Computes the Space Residual by assembling and solving the linear system.");
   options.addConfigOption< std::string >("SetNodalStatesCom","SetNodalStates Command to set the nodal states.");
   options.addConfigOption< std::string >("PrepareCom","Prepare Command to prepare the computation at each iteration. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetup Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
   options.addConfigOption< std::string >("SetupCom","Setup Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("ComputeTimeResidual","Computes the Time Residual by assembling and solving the linear system. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("BcComds","Types of the boundary conditions commands.");
   options.addConfigOption< std::vector<std::string> >("BcNames","Names for the configuration of the boundary conditions commands.");
}

//////////////////////////////////////////////////////////////////////////////

FiniteElementMethod::FiniteElementMethod(const std::string& name)
  : SpaceMethod(name),
    m_setup(),
    m_unSetup(),
    m_prepare(),
    m_extrapolateStates(),
    m_computeSpaceResidual(),
    m_computeTimeResidual(),
    m_inits(0),
    m_bcs(0)
{
  addConfigOptionsTo(this);
  m_data.reset(new FiniteElementMethodData(this));
  cf_assert(m_data.isNotNull());

  // set default value of builder for FiniteElementMethod
  // to be FEM_MeshDataBuilder
  m_builder = "FiniteElement";
  // set default global jacobian sparsity
  m_sparsity = "CellVertex";

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_prepareStr = "StdPrepare";
  setParameter("PrepareCom",&m_prepareStr);

  m_extrapolateStatesStr = "Null";
  setParameter("SetNodalStatesCom",&m_extrapolateStatesStr);

  m_computeSpaceResidualStr = "ExplicitComputeSpaceResCom";
  setParameter("ComputeSpaceResidual",&m_computeSpaceResidualStr);

  m_computeTimeResidualStr = "StdComputeTimeResCom";
  setParameter("ComputeTimeResidual",&m_computeTimeResidualStr);

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

FiniteElementMethod::~FiniteElementMethod()
{
  clearInitComs();
  clearBCComs();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::setCollaborator(MultiMethodHandle<LinearSystemSolver> lss)
{
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::clearInitComs()
{
  m_inits.clear();
  vector<SelfRegistPtr<FiniteElementMethodCom> >().swap(m_inits);
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::clearBCComs()
{
  m_bcs.clear();
  vector<SelfRegistPtr<FiniteElementMethodCom> >().swap(m_bcs);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> FiniteElementMethod::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<SpaceMethodData> FiniteElementMethod::getSpaceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::configure ( Config::ConfigArgs& args )
{
  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add here configures to the FiniteElementMethod

  configureCommand<FiniteElementMethodData,
                   FiniteElementMethodComProvider>
                       (args,
                        m_setup,
                        m_setupStr,
                        m_data);
  cf_assert(m_setup.isNotNull());

  configureCommand<FiniteElementMethodData,
                   FiniteElementMethodComProvider>
                       (args,
                        m_unSetup,
                        m_unSetupStr,
                        m_data);
  cf_assert(m_unSetup.isNotNull());

  configureCommand<FiniteElementMethodData,
                   FiniteElementMethodComProvider>
                       (args,
                        m_prepare,
                        m_prepareStr,
                        m_data);
  cf_assert(m_prepare.isNotNull());

  configureCommand<FiniteElementMethodData,
                   FiniteElementMethodComProvider>
    (args,
     m_extrapolateStates,
     m_extrapolateStatesStr,
     m_data);
  cf_assert(m_extrapolateStates.isNotNull());

  configureCommand<FiniteElementMethodData,
                   FiniteElementMethodComProvider>
                       (args,
                        m_computeSpaceResidual,
                        m_computeSpaceResidualStr,
                        m_data);
  cf_assert(m_computeSpaceResidual.isNotNull());


  configureCommand<FiniteElementMethodData,
                   FiniteElementMethodComProvider>
                       (args,
                        m_computeTimeResidual,
                        m_computeTimeResidualStr,
                        m_data);
  cf_assert(m_computeTimeResidual.isNotNull());

  clearInitComs();
  clearBCComs();

  cf_assert(m_initTypeStr.size() == m_initNameStr.size());

  m_inits.resize(m_initTypeStr.size());
  for(CFuint i = 0; i < m_inits.size(); ++i) {

    CFLog(INFO, "INIT type = " << m_initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << m_initNameStr[i] << "\n");

    configureCommand<FiniteElementMethodCom,
      FiniteElementMethodData,
      FiniteElementMethodComProvider>(args,m_inits[i],
              m_initTypeStr[i],
              m_initNameStr[i],
              m_data);
    cf_assert(m_inits[i].isNotNull());
  }

  cf_assert(m_bcTypeStr.size() == m_bcNameStr.size());

  m_bcs.resize(m_bcTypeStr.size());
  for(CFuint i = 0; i < m_bcs.size(); ++i) {

    CFLog(INFO, "BC type = " << m_bcTypeStr[i] << "\n");
    CFLog(INFO, "BC name = " << m_bcNameStr[i] << "\n");

    configureCommand<FiniteElementMethodCom,
      FiniteElementMethodData,
      FiniteElementMethodComProvider>(args,
              m_bcs[i],
              m_bcTypeStr[i],
              m_bcNameStr[i],
              m_data);
    cf_assert(m_bcs[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::initializeSolutionImpl(bool isRestart)
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

void FiniteElementMethod::prepareComputationImpl()
{
  CFAUTOTRACE;

  cf_assert(m_prepare.isNotNull());
  m_prepare->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;

  cf_assert(m_extrapolateStates.isNotNull());
  m_extrapolateStates->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(m_computeSpaceResidual.isNotNull());

  // set the residual factor for which residual and jacobian have to
  // be multiplied
  m_data->setResFactor(factor);

  m_computeSpaceResidual->execute();
  applyBC();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(m_computeTimeResidual.isNotNull());

  // set the residual factor for which residual and jacobian have to
  // be multiplied
  m_data->setResFactor(factor);

  m_computeTimeResidual->execute();
  checkMatrixFrozen();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::applyBCImpl()
{
  CFAUTOTRACE;
  m_data->setDirichletBCApplied(false);

  for(CFuint i = 0; i < m_bcs.size(); ++i)
  {
    cf_assert(m_bcs[i].isNotNull());
    m_bcs[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::setMethodImpl()
{
  CFAUTOTRACE;

  SpaceMethod::setMethodImpl();

  cf_assert(m_setup.isNotNull());
  m_setup->execute();

  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::unsetMethodImpl()
{
  CFAUTOTRACE;

  cf_assert(m_unSetup.isNotNull());
  m_unSetup->execute();

  unsetupCommandsAndStrategies();

  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::NumericalStrategy> > FiniteElementMethod::getStrategyList() const
{
  vector<Common::SafePtr<Framework::NumericalStrategy> > result;

  // add strategies here
  result.push_back(m_data->getResidualStrategy().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getJacobianStrategy().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getConvectiveTermComputer().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getDiffusiveTermComputer().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getIndepSourceTermComputer().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getLinearSourceTermComputer().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getInertiaTermComputer().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getConvectiveEntity().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getDiffusiveEntity().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getInertiaEntity().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getIndepSourceEntity().d_castTo<NumericalStrategy>());
  result.push_back(m_data->getLinearSourceEntity().d_castTo<NumericalStrategy>());

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethod::checkMatrixFrozen() const
{
  if(!m_data->isSysMatrixFrozen() && m_data->isSysMatrixFreezedEveryIteration()) {
    CFLog(INFO,"Freezing System Matrix" << "\n");
    m_data->setSysMatrixFrozen(true);
  }
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t FiniteElementMethod::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"FiniteElementMethod::beforeMeshUpdateAction() is DOING nothing!" << "\n");

  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t FiniteElementMethod::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"FiniteElementMethod::afterMeshUpdateAction() is DOING nothing!" << "\n");

  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
