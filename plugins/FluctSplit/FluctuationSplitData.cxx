#include "Common/CFPrintContainer.hh"

#include "Framework/ContourIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/OptionMethodStrategy.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"
#include "FluctSplit/ComputeJacobStrategy.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"
#include "FluctSplit/ComputeSourceTermFSM.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<FluctuationSplitData>,
          FluctuationSplitData, FluctSplitModule>
aNullFluctuationSplitComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOptionByType< OptionMethodStrategy< FluctuationSplitData, Splitter > >
     ("SysSplitter","System Splitter to distribute the residual.");

   options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature.");
   options.addConfigOption< std::string >("ContourIntegratorOrder","Order of the Integration to be used for numerical contour quadrature.");
   options.addConfigOption< std::string >("VolumeIntegratorOrder","Order of the Integration to be used for numerical volume quadrature.");
   options.addConfigOption< std::string >("FluctSplitStrategy","fluctuation splitting strategy");
   options.addConfigOption< std::string >("ArtDiffStrategy","Artificial diffusion strategy");
   options.addConfigOption< std::string >("JacobianStrategy","Jacobian computation strategy");
   options.addConfigOption< std::string >("ScalarSplitter","Scalar Splitter to distribute the residual.");
   options.addConfigOption< std::string >("JacobianSystemSplitter","System Splitter used to compute the Jacobian");
   options.addConfigOption< std::string >("JacobianScalarSplitter","Scalar Splitter used to compute the Jacobian");
   options.addConfigOption< std::string >("DistribVar","VarSet corresponding to Distribution variables.");
   options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration.");
   options.addConfigOption< std::string >("LinearVar","VarSet corresponding to Linearizing variables.");
   options.addConfigOption< std::string >("DiffusiveTerm","Diffusive term computer.");
   options.addConfigOption< vector<std::string> >("SourceTerm","Names of the source terms.");
   options.addConfigOption< vector<std::string> >("SourceTermSplitter","Type of scheme to distribute source terms.");
   options.addConfigOption< std::string >("JacobianFix", "Fix for the jacobian (against entropy violation or carbuncle).");
   options.addConfigOption< bool > ("isAxisymm","Tells if the SubSystem is axisymmetric.");
   options.addConfigOption< bool > ("includeSourceInFlux","Tells if a source term has to be added to the convective flux.");
   options.addConfigOption< bool > ("hasArtificialDiff","Tells if an artificial diffusion term should be added.");
   options.addConfigOption< bool>
     ("ScalarFirst","Flag telling if the scalar part has to be treated before the system part.");
}

//////////////////////////////////////////////////////////////////////////////

FluctuationSplitData::FluctuationSplitData(Common::SafePtr<Framework::Method> owner)
  : SpaceMethodData(owner),
    m_multipleSplitter(false),
    m_isOnlySysSplitter(false),
    m_rhsSysSplitter(),
    m_rhsSclSplitter(),
    m_jacobSysSplitter(),
    m_jacobSclSplitter(),
    m_sourceTermSplitter(),
    m_distribVar(),
    m_linearizer(),
    m_stComputer(),
    m_diffTermComputer(),
    m_jacobFixComputer(),
    m_solutionToDistMatTrans(),
    m_distToSolutionMatTrans(),
    m_linearToDistMatTrans(),
    m_solutionToLinearInUpdateMatTrans(),
    m_solutionToLinearMatTrans(),
    m_updateToLinearVecTrans(),
    m_updateToSolutionVecTrans(),
    m_solToUpdateInUpdateMatTrans(),
    m_updateToSolutionInUpdateMatTrans(),
    m_stdTrsGeoBuilder(),
    m_convergenceMtd(),
    m_linearVar(), // AL: possible memory leak: problems if you put this after m_distribVar
    m_distData(),
    m_resFactor(1.0),
    m_isInitializationPhase(false)
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  setParameter("SysSplitter", &m_rhsSysSplitter, std::string("Null"));
  getOptionT< OptionMethodStrategy< FluctuationSplitData, Splitter > >("SysSplitter")->putMethodData(this);

  m_scalarSplitterStr = "Null";
  setParameter("ScalarSplitter",&m_scalarSplitterStr) ;

  m_jacobSysSplitterStr = "Null";
  setParameter("JacobianSystemSplitter",&m_jacobSysSplitterStr);

  m_jacobSclSplitterStr = "Null";
  setParameter("JacobianScalarSplitter",&m_jacobSclSplitterStr);

  m_sourceTermSplitterStr = vector<std::string>();
  setParameter("SourceTermSplitter",&m_sourceTermSplitterStr);

  m_distribVarStr = "Null";
  setParameter("DistribVar",&m_distribVarStr);

  m_linearVarStr = "Null";
  setParameter("LinearVar",&m_linearVarStr);

  m_stComputerStr = vector<std::string>();
  setParameter("SourceTerm",&m_stComputerStr);

  m_diffTermComputerStr = "Null";
  setParameter("DiffusiveTerm",&m_diffTermComputerStr);

  m_jacobFixComputerStr = "Null";
  setParameter("JacobianFix",&m_jacobFixComputerStr);

  m_linearizerStr = "Null";
  m_solutionToDistMatTransStr = "Null";
  m_distToSolutionMatTransStr = "Null";
  m_linearToDistMatTransStr = "Null";
  m_solutionToLinearInUpdateMatTransStr = "Null";
  m_solutionToLinearMatTransStr = "Null";
  m_updateToLinearVecTransStr = "Null";
  m_updateToSolutionVecTransStr = "Null";
  m_solToUpdateInUpdateMatTransStr = "Null";
  m_updateToSolutionInUpdateMatTransStr = "Null";

  m_integratorOrderStr = "P1";
  setParameter("IntegratorOrder",&m_integratorOrderStr);

  m_cIntegratorOrderStr = "CFPolyOrder::MAXORDER";
  setParameter("ContourIntegratorOrder",&m_cIntegratorOrderStr);

  m_vIntegratorOrderStr = "CFPolyOrder::MAXORDER";
  setParameter("VolumeIntegratorOrder",&m_vIntegratorOrderStr);

  m_integratorQuadratureStr = "GaussLegendre";
  setParameter("IntegratorQuadrature",&m_integratorQuadratureStr);

  m_isAxisymm = false;
  setParameter("isAxisymm",&m_isAxisymm);

  m_includeSourceInFlux = false;
  setParameter("includeSourceInFlux",&m_includeSourceInFlux);

  m_fsStrategyName = "RD";
  setParameter("FluctSplitStrategy",&m_fsStrategyName);

  m_adStrategyName = "Null";
  setParameter("ArtDiffStrategy",&m_adStrategyName);

  m_jacobStrategyName = "Null";
  setParameter("JacobianStrategy",&m_jacobStrategyName);

  m_artdiff = false;
  setParameter("hasArtificialDiff",&m_artdiff); 
  
  m_scalarFirst = false;
  setParameter("ScalarFirst",&m_scalarFirst);
}

//////////////////////////////////////////////////////////////////////////////

FluctuationSplitData::~FluctuationSplitData()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"Configuring FluctuationSplitData\n");

  SpaceMethodData::configure(args);

  CFLog(VERBOSE,"Configuring VarSets\n");
  configureVarSets(args);

  CFLog(VERBOSE,"Configuring Splitters\n");
  configureSplitters(args);

  CFLog(VERBOSE,"Configuring SourceTermSplitters\n");
  configureSourceTermSplitters(args);

  CFLog(VERBOSE,"Configuring Linearizer\n");
  configureLinearizer(args);

  CFLog(VERBOSE,"Configuring FluctSplit Strategy\n");
  configureFluctSplitStrategy(args);

 CFLog(VERBOSE,"Configuring ArtDiff Strategy\n");
  configureArtDiffStrategy(args);

  CFLog(VERBOSE,"Configuring Jacobian Computation Strategy\n");
  configureJacobStrategy(args);

  // source term has to be configured
  // after the distribution variables have been configured
  CFLog(VERBOSE,"Configuring SourceTerm Computer\n");
  configureSourceTermComputer(args);

  CFLog(VERBOSE,"Configuring DiffusiveTerm Computer\n");
  configureDiffusiveTermComputer(args);

  CFLog(VERBOSE,"Configuring Transformers\n");
  configureTransformers(args);

  CFLog(VERBOSE,"Configuring Integrators\n");
  configureIntegrators(args);

  CFLog(INFO,"FluctuationSplit: Using Update VarSet: " << _updateVarStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Solution VarSet: " << _solutionVarStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Diffusive VarSet: " << _diffusiveVarStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Distribution VarSet: " << m_distribVarStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Linearizer: " << m_linearizerStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using ScalarSplitter: " << m_scalarSplitterStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using SystemSplitter: " << getSysSplitter()->getName() << "\n");
  CFLog(INFO,"FluctuationSplit: Using JacobianScalarSplitter: " << m_jacobSclSplitterStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using JacobianSystemSplitter: " << m_jacobSysSplitterStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using FluctuationSplitStrategy: " << getFluctSplitStrategyName() << "\n");
  CFLog(INFO,"FluctuationSplit: Using ArtificialDisffusionStrategy: " << getArtDiffStrategyName() << "\n");
  CFLog(INFO,"FluctuationSplit: Using MatrixTransformer Solution to Distribution Vars: " << m_solutionToDistMatTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using MatrixTransformer Distribution to Solution Vars: " << m_distToSolutionMatTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using MatrixTransformer Linear to Distribution Vars: " << m_linearToDistMatTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using MatrixTransformer Solution to Linear In Update Vars: " << m_solutionToLinearInUpdateMatTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using MatrixTransformer Solution to Linear Vars: " << m_solutionToLinearMatTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Transformer Update to Linear Vars: " << m_updateToLinearVecTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Transformer Update to Solution Vars: " << m_updateToSolutionVecTransStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Integrator Quadrature: " << m_integratorQuadratureStr << "\n");
  CFLog(INFO,"FluctuationSplit: Using Integrator Order: " << m_integratorOrderStr << "\n");
  CFLog(INFO, CFPrintContainer<vector<std::string> > ("FluctuationSplit: Using ComputeSourceTerm: ", &m_stComputerStr, m_stComputerStr.size()) << "\n");
  CFLog(INFO, CFPrintContainer<vector<std::string> > ("FluctuationSplit: Using SourceTermSplitter: ", &m_sourceTermSplitterStr, m_sourceTermSplitterStr.size()) << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureTransformers( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  Common::SafePtr<MatrixTransformer::PROVIDER> matTransProv = CFNULL;

  //--------------------------------------------------------------------------------------------
  // matrix transformer from solution to distribution vars
  m_solutionToDistMatTransStr = VarSetMatrixTransformer::getProviderName (physModel->getConvectiveName(), _solutionVarStr, m_distribVarStr, "Ref");

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_solutionToDistMatTransStr) )
    m_solutionToDistMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_solutionToDistMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_solutionToDistMatTrans = matTransProv->create(physModel->getImplementor());
  cf_assert(m_solutionToDistMatTrans.isNotNull());

  //--------------------------------------------------------------------------------------------
  // matrix transformer from distribution to solution vars
  matTransProv = CFNULL;
  m_distToSolutionMatTransStr = VarSetMatrixTransformer::getProviderName (physModel->getConvectiveName(), m_distribVarStr, _solutionVarStr, "Ref");

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_distToSolutionMatTransStr) )
    m_distToSolutionMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_distToSolutionMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_distToSolutionMatTrans = matTransProv->create(physModel->getImplementor());
  cf_assert(m_distToSolutionMatTrans.isNotNull());

  //--------------------------------------------------------------------------------------------
  // matrix transformer from linear to distribution vars
  matTransProv = CFNULL;
  m_linearToDistMatTransStr = VarSetMatrixTransformer::getProviderName (physModel->getConvectiveName(), m_linearVarStr, m_distribVarStr, "Ref");

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_linearToDistMatTransStr) )
    m_linearToDistMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_linearToDistMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_linearToDistMatTrans = matTransProv->create(physModel->getImplementor());
  cf_assert(m_linearToDistMatTrans.isNotNull());

  //--------------------------------------------------------------------------------------------
  // matrix transformer from linear to solution in solution vars
  // Usefull for instance in Picard Jacobian computation
  matTransProv = CFNULL;
  m_solutionToLinearInUpdateMatTransStr =  VarSetMatrixTransformer::getProviderName  (physModel->getConvectiveName(), _solutionVarStr,m_linearVarStr,_updateVarStr);

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_solutionToLinearInUpdateMatTransStr) )
    m_solutionToLinearInUpdateMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_solutionToLinearInUpdateMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_solutionToLinearInUpdateMatTrans = matTransProv->create(physModel->getImplementor());
  cf_assert(m_solutionToLinearInUpdateMatTrans.isNotNull());

  //--------------------------------------------------------------------------------------------
  // matrix transformer from solution to linear vars
  matTransProv = CFNULL;
  m_solutionToLinearMatTransStr = VarSetMatrixTransformer::getProviderName (physModel->getConvectiveName(), _solutionVarStr, m_linearVarStr, "Ref");

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_solutionToLinearMatTransStr) )
    m_solutionToLinearMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_solutionToLinearMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_solutionToLinearMatTrans = matTransProv->create(physModel->getImplementor());
  cf_assert(m_solutionToLinearMatTrans.isNotNull());

  //--------------------------------------------------------------------------------------------
  // matrix transformer from solution to update in update vars
  matTransProv = CFNULL;
  m_solToUpdateInUpdateMatTransStr = VarSetMatrixTransformer::getProviderName (physModel->getConvectiveName(), _solutionVarStr, _updateVarStr, _updateVarStr);

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_solToUpdateInUpdateMatTransStr) )
    m_solToUpdateInUpdateMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_solToUpdateInUpdateMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_solToUpdateInUpdateMatTrans.reset(matTransProv->create(physModel->getImplementor()));
  cf_assert(m_solToUpdateInUpdateMatTrans.getPtr() != CFNULL);

  //--------------------------------------------------------------------------------------------
  // matrix transformer from update to solution in update vars
  matTransProv = CFNULL;
  m_updateToSolutionInUpdateMatTransStr = VarSetMatrixTransformer::getProviderName (physModel->getConvectiveName(), _updateVarStr, _solutionVarStr, _updateVarStr);

  if ( ! Environment::Factory<MatrixTransformer>::getInstance().exists (m_updateToSolutionInUpdateMatTransStr) )
    m_updateToSolutionInUpdateMatTransStr = "Identity";

  matTransProv = Environment::Factory<MatrixTransformer>::getInstance().getProvider(m_updateToSolutionInUpdateMatTransStr);
  cf_assert(matTransProv.isNotNull());

  m_updateToSolutionInUpdateMatTrans.reset(matTransProv->create(physModel->getImplementor()));
  cf_assert(m_updateToSolutionInUpdateMatTrans.getPtr() != CFNULL);

  //--------------------------------------------------------------------------------------------
  // vector transformer from update to linear vars
  Common::SafePtr<VectorTransformer::PROVIDER> vecTransProv = CFNULL;
  m_updateToLinearVecTransStr = VarSetTransformer::getProviderName (physModel->getConvectiveName(), _updateVarStr,m_linearVarStr);

  if ( ! Environment::Factory<VectorTransformer>::getInstance().exists (m_updateToLinearVecTransStr) )
    m_updateToLinearVecTransStr = "Identity";

  vecTransProv = Environment::Factory<VectorTransformer>::getInstance().getProvider(m_updateToLinearVecTransStr);
  cf_assert(vecTransProv.isNotNull());

  m_updateToLinearVecTrans = vecTransProv->create(physModel->getImplementor());
  cf_assert(m_updateToLinearVecTrans.isNotNull());

  //--------------------------------------------------------------------------------------------
  // vector transformer from update to solution vars
  vecTransProv = CFNULL;
  m_updateToSolutionVecTransStr = VarSetTransformer::getProviderName (physModel->getConvectiveName(), _updateVarStr, _solutionVarStr);

  if ( ! Environment::Factory<VectorTransformer>::getInstance().exists (m_updateToSolutionVecTransStr) )
    m_updateToSolutionVecTransStr = "Identity";

  vecTransProv = Environment::Factory<VectorTransformer>::getInstance().getProvider(m_updateToSolutionVecTransStr);
  cf_assert(vecTransProv.isNotNull());

  m_updateToSolutionVecTrans = vecTransProv->create(physModel->getImplementor());
  cf_assert(m_updateToSolutionVecTrans.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureSplitters ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

//   CFLogDebugMed("Configuring System Splitter: " << m_sysSplitterStr << "\n");
//
//   SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,Splitter> > sys_prov =
//     Environment::Factory<Splitter>::getInstance().getProvider(m_sysSplitterStr);
//   cf_assert(sys_prov.isNotNull());
//
//   m_rhsSysSplitter = sys_prov->create(m_sysSplitterStr, SharedPtr<FluctuationSplitData>(this));
//   configureNested ( m_rhsSysSplitter.getPtr(), args );
//   cf_assert(m_rhsSysSplitter.isNotNull());

  configureNested ( m_rhsSysSplitter.getPtr(), args );

  CFLogDebugMed("Configuring Scalar Splitter: " << m_scalarSplitterStr << "\n");
  SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,Splitter> > scal_prov =
    Environment::Factory<Splitter>::getInstance().getProvider(m_scalarSplitterStr);
  cf_assert(scal_prov.isNotNull());

  m_rhsSclSplitter = scal_prov->create(m_scalarSplitterStr, SharedPtr<FluctuationSplitData>(this));
  configureNested ( m_rhsSclSplitter.getPtr(), args );
  cf_assert(m_rhsSclSplitter.isNotNull());

  // configure a separate jacobian splitter only if user supplied it
  // if not, use the same one as the rhs
  if (m_jacobSysSplitterStr == "Null")
  {
    m_jacobSysSplitter = m_rhsSysSplitter;
  }
  else
  {
    CFLogDebugMed("Configuring Jacobian System Splitter: " << m_jacobSysSplitterStr << "\n");
    SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,Splitter> > jacob_prov =
      Environment::Factory<Splitter>::getInstance().getProvider(m_jacobSysSplitterStr);
    cf_assert(jacob_prov.isNotNull());

    m_jacobSysSplitter = jacob_prov->create(m_jacobSysSplitterStr, SharedPtr<FluctuationSplitData>(this));
    configureNested ( m_jacobSysSplitter.getPtr(), args );
    cf_assert(m_jacobSysSplitter.isNotNull());
  }

  // configure a separate jacobian splitter only if user supplied it
  // if not, use the same one as the rhs
  if (m_jacobSclSplitterStr == "Null")
  {
    m_jacobSclSplitter = m_rhsSclSplitter;
  }
  else
  {
    CFLogDebugMed("Configuring Jacobian Scalar Splitter: " << m_jacobSclSplitterStr << "\n");
    SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,Splitter> > jacob_scal_prov =
      Environment::Factory<Splitter>::getInstance().getProvider(m_jacobSclSplitterStr);
    cf_assert(jacob_scal_prov.isNotNull());

    m_jacobSclSplitter = jacob_scal_prov->create(m_jacobSclSplitterStr, SharedPtr<FluctuationSplitData>(this));
    configureNested ( m_jacobSclSplitter.getPtr(), args );
    cf_assert(m_jacobSclSplitter.isNotNull());
  }


  // setting the flags for system/scalar splitting
  m_multipleSplitter = (!(m_rhsSysSplitter->isNull()) &&  !(m_rhsSclSplitter->isNull()));
  m_isOnlySysSplitter =  (m_rhsSclSplitter->isNull()) && (!(m_rhsSysSplitter->isNull()));

  // by default the splitters are set to the
  // explicit computation
  m_SysSplitter = m_rhsSysSplitter.getPtr();
  m_SclSplitter = m_rhsSclSplitter.getPtr();

  // jacobian fix computer
  SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,ComputeJacobianFix> > jfprov =
    Environment::Factory<ComputeJacobianFix>::getInstance().getProvider(m_jacobFixComputerStr);
  cf_assert(jfprov.isNotNull());

  m_jacobFixComputer = jfprov->create(m_jacobFixComputerStr, SharedPtr<FluctuationSplitData>(this));
  cf_assert(m_jacobFixComputer.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureSourceTermSplitters ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFLogDebugMed("Configuring SourceTermSplitters\n");

  if (m_sourceTermSplitterStr.size() == 0) {
    m_sourceTermSplitterStr.resize(1);
    m_sourceTermSplitterStr[0] = "Null";
  }
  cf_assert(m_sourceTermSplitterStr.size() > 0);
  m_sourceTermSplitter.resize(m_sourceTermSplitterStr.size());
  cf_assert(m_sourceTermSplitter.size() == m_sourceTermSplitterStr.size());

  for (CFuint i = 0; i < m_sourceTermSplitterStr.size(); ++i) {
    SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,SourceTermSplitter> > provST =
      Environment::Factory<SourceTermSplitter>::getInstance().getProvider(m_sourceTermSplitterStr[i]);
    cf_assert(provST.isNotNull());

    m_sourceTermSplitter[i] = provST->create
      (m_sourceTermSplitterStr[i], SharedPtr<FluctuationSplitData>(this));
    cf_assert(m_sourceTermSplitter[i].isNotNull());
    configureNested(m_sourceTermSplitter[i].getPtr(), args);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureLinearizer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  m_linearizerStr = physModel->getConvectiveName()
                                    + "Linear"
                                    + m_linearVarStr;

  CFLogDebugMed("Configuring JacobianLinearizer: " << m_linearizerStr << "\n");

  Common::SafePtr<JacobianLinearizer::PROVIDER> lin_prov =
    Environment::Factory<JacobianLinearizer>::getInstance().getProvider(m_linearizerStr);

  cf_assert(lin_prov.isNotNull());

  m_linearizer = lin_prov->create(physModel);

  cf_assert(m_linearizer.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureFluctSplitStrategy( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"FluctuationSplitStrategy : " << getFluctSplitStrategyName() << "\n");

  SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,FluctuationSplitStrategy> > prov;
  try
  {
    prov = Environment::Factory<FluctuationSplitStrategy>::getInstance().
           getProvider(getFluctSplitStrategyName());
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    throw;
  }

  cf_assert(prov.isNotNull());
  m_fsStrategy = prov->create(getFluctSplitStrategyName(),SharedPtr<FluctuationSplitData>(this));
  cf_assert(m_fsStrategy.isNotNull());
  configureNested ( m_fsStrategy.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureArtDiffStrategy( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  configureStrategy < ArtificialDiffusionStrategy,
                      FluctuationSplitData >
      ( args, m_adStrategy, getArtDiffStrategyName(), getArtDiffStrategyName(), SharedPtr<FluctuationSplitData>(this));
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureJacobStrategy( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"ComputeJacobStrategy : " << getComputeJacobianStrategyName() << "\n");

  SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,ComputeJacobStrategy> > prov = CFNULL;
  try
  {
    prov = Environment::Factory<ComputeJacobStrategy>::getInstance().
           getProvider(getComputeJacobianStrategyName());
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    throw;
  }

  cf_assert(prov.isNotNull());
  m_jacobStrategy = prov->create(getComputeJacobianStrategyName(),SharedPtr<FluctuationSplitData>(this));
  cf_assert(m_jacobStrategy.isNotNull());
  configureNested ( m_jacobStrategy.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureSourceTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  CFLogDebugMed("Configuring SourceTerm Computer\n");

  if (m_stComputerStr.size() == 0) {
    m_stComputerStr.resize(1);
    m_stComputerStr[0] = "Null";
  }


  cf_assert(m_stComputerStr.size() > 0);
  m_stComputer.resize(m_sourceTermSplitterStr.size());
  cf_assert(m_stComputer.size() == m_stComputerStr.size());

  SafePtr<FluctuationSplitData> thisPtr(this);

  for (CFuint i = 0; i < m_stComputerStr.size(); ++i) {
    const std::string stName = m_stComputerStr[i];
    m_stComputer[i] = Environment::Factory<ComputeSourceTermFSM>::getInstance().
      getProvider(stName)->create(stName, thisPtr);
    configureNested(m_stComputer[i].getPtr(),  args );
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureVarSets( Config::ConfigArgs& args )
 {
  CFAUTOTRACE;

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  CFLogDebugMed("Configuring Distribute VarSet\n");
  CFLogDebugMed("VarSet is : " << physModel->getConvectiveName() + m_distribVarStr << "\n");

  m_distribVar.reset(Environment::Factory<ConvectiveVarSet>::getInstance().
    getProvider(physModel->getConvectiveName() + m_distribVarStr)->create(physModel->getImplementor()->getConvectiveTerm()));

  cf_assert(m_distribVar.isNotNull());

  CFLogDebugMed("Configuring Linearization VarSet\n");
  CFLogDebugMed("VarSet is : " << physModel->getConvectiveName() + m_linearVarStr << "\n");

  m_linearVar.reset(Environment::Factory<ConvectiveVarSet>::getInstance().
        getProvider(physModel->getConvectiveName() + m_linearVarStr)->
        create(physModel->getImplementor()->getConvectiveTerm()));

  cf_assert(m_linearVar.isNotNull());
 }

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureIntegrators( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFPolyOrder::Type order;

  CFLog(VERBOSE,"Configuring ContourIntegrator\n");

  const CFQuadrature::Type quadType =
    CFQuadrature::Convert::to_enum( m_integratorQuadratureStr );

  const CFPolyOrder::Type orderGen =
    CFPolyOrder::Convert::to_enum( m_integratorOrderStr );

  // either use the specific for the contour or the general order, which is also the default
  order = CFPolyOrder::Convert::to_enum( m_cIntegratorOrderStr );
  if ( order == CFPolyOrder::MAXORDER ) { order = orderGen; }

  try {
     m_ContourIntegrator.setIntegrationForAllGeo(quadType,order);
  }
  catch (NoSuchValueException& e) {
    // output what went wrong
    CFLog(ERROR, e.what() << "\n");
    CFLog(ERROR, "Integrator Quadrature : " << m_integratorQuadratureStr << "\n");
    CFLog(ERROR, "Integrator Order : " << m_integratorOrderStr << "\n");
    throw; // rethrow exception
  }

  CFLog(VERBOSE,"Configuring VolumeIntegrator\n");

  // either use the specific for the contour or the general order, which is also the default
  order = CFPolyOrder::Convert::to_enum( m_vIntegratorOrderStr );
  if ( order == CFPolyOrder::MAXORDER ) { order = orderGen; }

  try {
     m_VolumeIntegrator.setIntegrationForAllGeo(quadType,order);
  }
  catch (NoSuchValueException& e) {
    // output what went wrong
    CFLog(ERROR, e.what() << "\n");
    CFLog(ERROR, "Integrator Quadrature : " << m_integratorQuadratureStr << "\n");
    CFLog(ERROR, "Integrator Order : " << m_integratorOrderStr << "\n");
    throw; // rethrow exception
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::configureDiffusiveTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SafePtr<BaseMethodStrategyProvider<FluctuationSplitData,ComputeDiffusiveTerm> > prov;
  try {
    prov = Environment::Factory<ComputeDiffusiveTerm>::getInstance().
           getProvider(m_diffTermComputerStr);
  }
  catch (Common::NoSuchValueException& e) {
    CFLog(VERBOSE, e.what() << "\n");
    throw;
  }

  CFLogDebugMed("Configuring ComputeDiffusiveTerm"<< m_diffTermComputerStr << "\n");

  cf_assert(prov.isNotNull());

  m_diffTermComputer = prov->create(m_diffTermComputerStr,SharedPtr<FluctuationSplitData>(this));
  cf_assert(m_diffTermComputer.isNotNull());
  configureNested ( m_diffTermComputer.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

CFuint FluctuationSplitData::getBlockSeparator() const
{

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  if (m_rhsSysSplitter->getName() != "Null" && m_scalarSplitterStr != "Null") {
    return m_distribVar->getBlockSeparator();
  }
  return physModel->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::setup()
{
  CFAUTOTRACE;

  SpaceMethodData::setup();

  // setup integrators
  m_ContourIntegrator.setup();
  m_VolumeIntegrator.setup();

  // this should be moved away => FluctuationSplit.cxx
  m_jacobSclSplitter->setup();
  m_jacobSysSplitter->setup();

  // setup TRS Geo builder
  m_stdTrsGeoBuilder.setup();

  getDistribVar()->setup();
  getUpdateVar()->setup();
  getLinearVar()->setup();
  getDiffusiveVar()->setup();
  getSolutionVar()->setup();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  getSolutionToDistribMatTrans()->setup(maxNbStatesInCell);
  getDistribToSolutionMatTrans()->setup(maxNbStatesInCell);
  getLinearToDistribMatTrans()->setup(maxNbStatesInCell);
  getSolutionToLinearInUpdateMatTrans()->setup(maxNbStatesInCell);
  getSolutionToLinearMatTrans()->setup(maxNbStatesInCell);
  getUpdateToLinearVecTrans()->setup(maxNbStatesInCell);
  getSolToUpdateInUpdateMatTrans()->setup(maxNbStatesInCell);
  getUpdateToSolutionInUpdateMatTrans()->setup(maxNbStatesInCell);
  getUpdateToSolutionVecTrans()->setup(maxNbStatesInCell);
  getDiffusiveTermComputer()->setMeshData();
  getDiffusiveTermComputer()->setDiffusiveVarSet(getDiffusiveVar());
  getDiffusiveTermComputer()->setUpdateVarSet(getUpdateVar());

  // set the maximum number of states in a cell
  // it will be useful to resize the LinearizationData in the
  // concrete physical model
  getLinearizer()->setMaxNbStates(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());

  // set up the distribution data
  getDistributionData().setup();

}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitData::unsetup()
{
  SpaceMethodData::unsetup();
  
  // unsetup TRS Geo builder
  m_stdTrsGeoBuilder.unsetup();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<FluctuationSplitStrategy> FluctuationSplitData::getFluctSplitStrategy() const
{
  cf_assert(m_fsStrategy.isNotNull());
  return m_fsStrategy.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ArtificialDiffusionStrategy> FluctuationSplitData::getArtificialDiffusionStrategy() const
{
  cf_assert(m_adStrategy.isNotNull());
  return m_adStrategy.getPtr();
}
//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ComputeJacobStrategy> FluctuationSplitData::getJacobStrategy() const
{
  cf_assert(m_jacobStrategy.isNotNull());
  return m_jacobStrategy.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<VolumeIntegrator> FluctuationSplitData::getVolumeIntegrator()
{
  return &m_VolumeIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ContourIntegrator> FluctuationSplitData::getContourIntegrator()
{
  return &m_ContourIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<vector<SelfRegistPtr<ComputeSourceTermFSM> > > FluctuationSplitData::
getSourceTermComputer()
{
  cf_assert(m_stComputer.size() > 0);
  return &m_stComputer;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<ComputeSourceTermFSM> FluctuationSplitData::getSourceTermComputer(CFuint is)
{
  cf_assert(is < m_stComputer.size());
  return m_stComputer[is].getPtr();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

