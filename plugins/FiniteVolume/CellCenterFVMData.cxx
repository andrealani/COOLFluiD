#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/ComputeDiffusiveFlux.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "FiniteVolume/FVMCC_PolyRec.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<CellCenterFVMData>, CellCenterFVMData, FiniteVolumeModule> nullCellCenterFVMComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Limiter","Limiter to be used.");
  options.addConfigOption< std::string >("PolyRec","PolyReconstructor to be used.");
  options.addConfigOption< std::string >("DiffusiveFlux","Diffusive flux computer.");
  options.addConfigOption< std::vector<std::string> >("EquationFilter","Equation filter type(s).");
  options.addConfigOption< std::vector<std::string> >("SourceTerm","Source term type(s).");
  options.addConfigOption< std::vector<std::string> >("TRSsWithGhostsOnFace","TRSs on which ghost states lay on the face."); 
  options.addConfigOption< std::vector<std::string> >("TRSsWithNoBC","TRSs for which no BC must be applied.");
  options.addConfigOption< std::string >("DerivativeStrategy","Derivative computation strategy.");
  options.addConfigOption< std::string >("LinearVar","Linearization variable set."); 
  options.addConfigOption< std::string >("ReconstrVar","Reconstruction variable set.");
  options.addConfigOption< std::string >("NodalExtrapolation","Nodal extrapolation strategy.");
  options.addConfigOption< bool >("isAxisymm","Tells if the simulation is axisymmetric.");
  options.addConfigOption< std::string >("GeoDataComputer","Computer of geometric data (e.g. normals, volumes).");
  options.addConfigOption< std::string >("FluxSplitter","FluxSplitter to compute flux.");
  options.addConfigOption< bool >("UseAnalyticalConvJacob","Use the analytical jacobian of convective fluxes");
  options.addConfigOption< bool >("ReconstructSolutionVars", "Reconstruct the solution variables instead of the update ones");
  options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature.");
  options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration.");

}
      
//////////////////////////////////////////////////////////////////////////////

CellCenterFVMData::CellCenterFVMData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  _isInitializationPhase(false),
  _isPreProcessedSolution(false),
  _volumeIntegrator(),
  _diffusiveFlux(),
  _geoDataComputer(),
  _fluxSplitter(),
  _polyRec(),
  _limiter(),
  _nStatesExtrapolator(),
  _convergenceMtd(),
  _eqFilters(),
  _stComputer(),
  _solToUpdateInUpdateMatTrans(),
  _updateToSolutionInUpdateMatTrans(),
  _updateToSolutionVecTrans(),
  _solutionToLinearVecTrans(), 
  _updateToReconstrVecTrans(),
  _reconstrToUpdateVecTrans(),
  _linearizer(),
  _faceTrsGeoBuilder(),
  _faceCellTrsGeoBuilder(),
  _cellTrsGeoBuilder(),
  _geoWithNodesBuilder(),
  _currFace(CFNULL),
  _bcMap(),
  _unitNormal(),
  _preProcessBCFlag(false),
  _useAverageFlux(false),
  _hasSourceTerm(false),
  _buildAllCells(false),
  _resFactor(1.0)
{
  addConfigOptionsTo(this);

  _updateVarStr = "Prim";
  _solutionVarStr = "Prim";
  
  _linearVarStr = "Null";
  setParameter("LinearVar",&_linearVarStr);
   
  _reconstructVarStr = "Null";
  setParameter("ReconstrVar",&_reconstructVarStr);
  
  _geoDataComputerStr = "FVMCC";
  setParameter("GeoDataComputer",&_geoDataComputerStr);
  
  _fluxSplitterStr = "Null";
  setParameter("FluxSplitter",&_fluxSplitterStr);
  
  _diffusiveFluxStr = "Null";
  setParameter("DiffusiveFlux",&_diffusiveFluxStr);
  
  setParameter("DerivativeStrategy",&_derivComputerStr);
  _derivComputerStr = "Null";
  
  _polyRecStr = "Constant";
  setParameter("PolyRec",&_polyRecStr);

  _limiterStr = "Null";
  setParameter("Limiter",&_limiterStr);
  
  _integratorOrderStr = "P1";
  setParameter("IntegratorOrder",&_integratorOrderStr);
  
  _integratorQuadratureStr = "GaussLegendre";
  setParameter("IntegratorQuadrature",&_integratorQuadratureStr);
  
  _nStatesExtrapolatorStr = "DistanceBased";
  setParameter("NodalExtrapolation",&_nStatesExtrapolatorStr);
  
  _eqFiltersStr = vector<std::string>();
  setParameter("EquationFilter",&_eqFiltersStr);
  
  _stComputerStr = vector<std::string>();
  setParameter("SourceTerm",&_stComputerStr);
  
  _trssWithGhostsOnFace = vector<std::string>();
  setParameter("TRSsWithGhostsOnFace",&_trssWithGhostsOnFace);

  _trssWithNoBC = vector<std::string>();
  setParameter("TRSsWithNoBC",&_trssWithNoBC);
  
  _isAxisymm = false;
  setParameter("isAxisymm",&_isAxisymm);

  _useAnalyticalConvJacob = false;
  setParameter("UseAnalyticalConvJacob",&_useAnalyticalConvJacob);

  _reconstructSolVars = false;
  setParameter("ReconstructSolutionVars",&_reconstructSolVars);
}

//////////////////////////////////////////////////////////////////////////////

CellCenterFVMData::~CellCenterFVMData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "CellCenterFVMData::configure() => before SpaceMethodData::configure()\n");
  
  SpaceMethodData::configure(args);
  
  // AL: volume integrator is still needed for ALE
  const CFQuadrature::Type quadType = CFQuadrature::Convert::to_enum(_integratorQuadratureStr);
  const CFPolyOrder::Type order = CFPolyOrder::Convert::to_enum( _integratorOrderStr );
  _volumeIntegrator.setFactoryRegistry(getFactoryRegistry());
  _volumeIntegrator.setIntegrationForAllGeo(quadType,order);
  
  CFLog(VERBOSE, "CellCenterFVMData::configure() => after SpaceMethodData::configure()\n");
  
  // default value for the linear variables is the update variables
  if (_linearVarStr == "Null") {
    _linearVarStr = _updateVarStr;
  }
  
  if (_reconstructVarStr == "Null") {
    _reconstructVarStr = _updateVarStr;
  }
  
  CFLog(INFO,"CellCenterFVM: configureGeoDataComputer()\n");
  configureGeoDataComputer(args);
  
  CFLog(INFO,"CellCenterFVM: configureFluxSplitter()\n");
  configureFluxSplitter(args);
  
  CFLog(INFO,"CellCenterFVM: configureDiffusiveFluxComputer()\n");
  configureDiffusiveFluxComputer(args);
  
  CFLog(INFO,"CellCenterFVM: configureVarSetTransformers()\n");
  configureVarSetTransformers(args);
  
  CFLog(INFO,"CellCenterFVM: configureJacobianLinearizer()\n");
  configureJacobianLinearizer(args);
  
  CFLog(INFO,"CellCenterFVM: configurePolyReconstructor()\n");
  configurePolyReconstructor(args);

  CFLog(INFO,"CellCenterFVM: configureLimiter()\n");
  configureLimiter(args);

  CFLog(INFO,"CellCenterFVM: configureNodalStatesExtrapolator()\n");
  configureNodalStatesExtrapolator(args);

  CFLog(INFO,"CellCenterFVM: configureEquationFilters()\n");
  configureEquationFilters(args);
  
  CFLog(INFO,"CellCenterFVM: configureSourceTermComputer()\n");
  configureSourceTermComputer(args);
  
  CFLog(INFO,"CellCenterFVM: Using GeoDataComputer: " << _geoDataComputerStr << "\n");
  CFLog(INFO,"CellCenterFVM: Using FluxSplitter: " << _fluxSplitterStr << "\n");
  CFLog(INFO,"CellCenterFVM: Using Update VarSet: " << _updateVarStr << "\n");
  CFLog(INFO,"CellCenterFVM: Using Solution VarSet: " << _solutionVarStr << "\n");
  CFLog(INFO,"CellCenterFVM: Using Diffusive VarSet: " << _diffusiveVarStr << "\n");
  CFLog(INFO,"CellCenterFVM: Using Linear VarSet: " << _linearVarStr << "\n");
  CFLog(INFO,"CellCenterFVM: Using NodalStatesExtrapolator: " <<  _nStatesExtrapolatorStr << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureGeoDataComputer ( Config::ConfigArgs& args )
{
  std::string name = _geoDataComputerStr;
  
  CFLogDebugMin("CellCenterFVM: Using GeoDataComputer: " << name << "\n");
  
  SharedPtr<CellCenterFVMData> thisPtr(this);
  Common::SafePtr<BaseMethodStrategyProvider<CellCenterFVMData,GeoDataComputer<CellCenterFVMData> > > prov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), GeoDataComputer<CellCenterFVMData>, name);
  
  cf_assert(prov.isNotNull());
  
  _geoDataComputer = prov->create(name,thisPtr);
  _geoDataComputer->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _geoDataComputer.getPtr(), args ); 
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureFluxSplitter ( Config::ConfigArgs& args )
{
  std::string name = _fluxSplitterStr;
  
  CFLogDebugMin("CellCenterFVM: Using FluxSplitter: " << name << "\n");

  SharedPtr<CellCenterFVMData> thisPtr(this);
  Common::SafePtr<BaseMethodStrategyProvider<CellCenterFVMData,FluxSplitter<CellCenterFVMData> > > prov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), FluxSplitter<CellCenterFVMData>, name);

  cf_assert(prov.isNotNull());

  _fluxSplitter = prov->create(name,thisPtr);
  _fluxSplitter->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _fluxSplitter.getPtr(), args ); 
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureDiffusiveFluxComputer ( Config::ConfigArgs& args )
{
  SharedPtr<CellCenterFVMData> thisPtr(this);

  Common::SafePtr<BaseMethodStrategyProvider<CellCenterFVMData,
     ComputeDiffusiveFlux > > prov = CFNULL;

  try {
    prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeDiffusiveFlux, _diffusiveFluxStr);
  }
  catch (Common::NoSuchValueException& e) {
   _diffusiveFluxStr = "Null";

    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing NullDiffusiveFlux instead ..." << "\n");
    prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeDiffusiveFlux, _diffusiveFluxStr);
  }

  cf_assert(prov.isNotNull());
  _diffusiveFlux = prov->create(_diffusiveFluxStr, thisPtr);
  _diffusiveFlux->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _diffusiveFlux.getPtr(), args );

  // configuration of the DerivativeComputer

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  if (_derivComputerStr == "Null" && _fluxSplitterStr != "Null") {
    if (physModel->getDim() >= 2) {
      const std::string dim = (physModel->getDim() == 2) ? "2D" : "3D";
      _derivComputerStr = "DiamondVolume" + dim;
    }
  }

  if (_diffusiveFluxStr == "Null") {_derivComputerStr = "Null";}

  CFLog(DEBUG_MED,"configuringDerivativeComputer()\n");
  CFLog(DEBUG_MED,"DerivativeComputer: " << _derivComputerStr << "\n");

  Common::SafePtr<BaseMethodStrategyProvider<CellCenterFVMData,
     DerivativeComputer > > dcprov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), DerivativeComputer, _derivComputerStr);
  cf_assert(dcprov.isNotNull());
  _derivComputer = dcprov->create(_derivComputerStr,thisPtr);
  
  configureNested ( _derivComputer.getPtr(), args );
  _derivComputer->setFactoryRegistry(getFactoryRegistry());
  cf_assert(_derivComputer.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configurePolyReconstructor ( Config::ConfigArgs& args )
{
  
  SharedPtr<CellCenterFVMData> thisPtr(this);

  Common::SafePtr<BaseMethodStrategyProvider<CellCenterFVMData,
    PolyReconstructor<CellCenterFVMData> > > prov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), PolyReconstructor<CellCenterFVMData>, _polyRecStr);
  cf_assert(prov.isNotNull());
  _polyRec = prov->create(_polyRecStr,thisPtr).d_castTo<FVMCC_PolyRec>();
  _polyRec->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _polyRec.getPtr(), args );
  
  cf_assert(_polyRec.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureLimiter ( Config::ConfigArgs& args )
{

  SharedPtr<CellCenterFVMData> thisPtr(this);

  Common::SafePtr<BaseMethodStrategyProvider<CellCenterFVMData,Limiter<CellCenterFVMData> > > prov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), Limiter<CellCenterFVMData>, _limiterStr);
  cf_assert(prov.isNotNull());
  _limiter = prov->create(_limiterStr,thisPtr);
  _limiter->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _limiter.getPtr(), args );
  
  cf_assert(_limiter.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureNodalStatesExtrapolator ( Config::ConfigArgs& args )
{
  // create the nodal states extrapolator
  SharedPtr<CellCenterFVMData> thisPtr(this);

  Common::SafePtr<BaseMethodStrategyProvider<
    CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData> > > prov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), NodalStatesExtrapolator<CellCenterFVMData>,
			 _nStatesExtrapolatorStr);
  
  cf_assert(prov.isNotNull());
  _nStatesExtrapolator = prov->create(_nStatesExtrapolatorStr,thisPtr);
  _nStatesExtrapolator->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _nStatesExtrapolator.getPtr(), args );

  cf_assert(_nStatesExtrapolator.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureEquationFilters ( Config::ConfigArgs& args )
{
  SharedPtr<CellCenterFVMData> thisPtr(this);
  
  if (_eqFiltersStr.size() > 0) {
    // resize the vector of equation filters
    _eqFilters.resize(_eqFiltersStr.size());
    
    Common::SafePtr<BaseMethodStrategyProvider
		    <CellCenterFVMData,EquationFilter<CellCenterFVMData> > > prov = CFNULL;
    
    for (CFuint i = 0; i < _eqFiltersStr.size(); ++i) {
      try {
	prov = 
	  FACTORY_GET_PROVIDER(getFactoryRegistry(), EquationFilter<CellCenterFVMData>, _eqFiltersStr[i]);
      }
      catch (Common::NoSuchValueException& e ) {
	CFLog(VERBOSE, e.what() << "\n");
	CFLog(VERBOSE, "Choosing NullEquationFilter instead ..." << "\n");
	_eqFiltersStr[i] = "Null";
	
	prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), EquationFilter<CellCenterFVMData>, _eqFiltersStr[i]);
      }
      
      _eqFilters[i] = prov->create(_eqFiltersStr[i], thisPtr);
    }
  }
  else {
    // create a single Null equatiion filter
    _eqFilters.resize(1);
    _eqFilters[0] = FACTORY_GET_PROVIDER(getFactoryRegistry(), EquationFilter<CellCenterFVMData>, "Null")->create("Null",thisPtr);
  }
  
  for (CFuint i = 0; i < _eqFilters.size(); ++i) {
    configureNested(_eqFilters[i].getPtr(), args);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureSourceTermComputer ( Config::ConfigArgs& args )
{
  SharedPtr<CellCenterFVMData> thisPtr(this);

  // this will prevent from even calling the object to compute the source term
  // if the source term is Null
  _hasSourceTerm = (_isAxisymm || _stComputerStr.size() > 0 ||
		    _solutionVar->hasSourceTerm()) ? true : false;

  if (_hasSourceTerm) {
    try {
      if (_stComputerStr.size() == 0) {
        throw NoSuchValueException
          (FromHere(), "CellCenterFVMData::configureSourceTermComputer() => term not found");
      }
    }
    catch (Common::NoSuchValueException& e ) {
      CFout << e.what() << "\n";
    }

    // resize the vector of source terms
    _stComputer.resize(_stComputerStr.size());
    
    Common::SafePtr<BaseMethodStrategyProvider
		    <CellCenterFVMData,
		     ComputeSourceTerm<CellCenterFVMData> > > prov = CFNULL;
    
    for (CFuint i = 0; i < _stComputerStr.size(); ++i) {
      try {
        prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeSourceTerm<CellCenterFVMData>, _stComputerStr[i]);
      }
      catch (Common::NoSuchValueException& e ) {
        CFLog(VERBOSE, e.what() << "\n");
        CFLog(VERBOSE, "Choosing NullComputeSourceTerm instead ..." << "\n");
        _stComputerStr[i] = "Null";

        prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeSourceTerm<CellCenterFVMData>, _stComputerStr[i]);
      }

      _stComputer[i] = prov->create(_stComputerStr[i], thisPtr);
    }
  }
  else {
    // create a single Null source term
    _stComputer.resize(1);
    _stComputer[0] = FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeSourceTerm<CellCenterFVMData>, "Null")->create("Null",thisPtr);
  }

  for (CFuint i = 0; i < _stComputerStr.size(); ++i) {
    _stComputer[i]->setFactoryRegistry(getFactoryRegistry());
    configureNested(_stComputer[i].getPtr(), args);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureVarSetTransformers ( Config::ConfigArgs& args )
{
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string solToUpdateInUpdateMatTransStr =
    VarSetMatrixTransformer::getProviderName
    (physModel->getConvectiveName(), _solutionVarStr,
     _updateVarStr, _updateVarStr);

  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
	solToUpdateInUpdateMatTransStr << "\n");

  Common::SafePtr<VarSetMatrixTransformer::PROVIDER> matTransProv = CFNULL;

  try {
    matTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetMatrixTransformer, solToUpdateInUpdateMatTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    solToUpdateInUpdateMatTransStr = "Identity";

    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetMatrixTransformer instead ..." << "\n");
    matTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetMatrixTransformer, solToUpdateInUpdateMatTransStr);
  }

  cf_assert(matTransProv.isNotNull());
  _solToUpdateInUpdateMatTrans.reset(matTransProv->create(physModel->getImplementor()));
  cf_assert(_solToUpdateInUpdateMatTrans.getPtr() != CFNULL);
    
  std::string updateToSolutionInUpdateMatTransStr =
    VarSetMatrixTransformer::getProviderName
    (physModel->getConvectiveName(),_updateVarStr,_solutionVarStr,_updateVarStr);

  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
	updateToSolutionInUpdateMatTransStr << "\n");
  try {
    matTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetMatrixTransformer, updateToSolutionInUpdateMatTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    updateToSolutionInUpdateMatTransStr = "Identity";
    
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetMatrixTransformer instead ..." << "\n");
    matTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetMatrixTransformer, updateToSolutionInUpdateMatTransStr);
  }
  
  cf_assert(matTransProv.isNotNull());
  _updateToSolutionInUpdateMatTrans.reset
    (matTransProv->create(physModel->getImplementor()));
  cf_assert(_updateToSolutionInUpdateMatTrans.getPtr() != CFNULL);
  
  Common::SafePtr<VarSetTransformer::PROVIDER> vecTransProv = CFNULL;
  
  std::string updateToSolutionVecTransStr =
    VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), _updateVarStr, _solutionVarStr);

  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
	updateToSolutionVecTransStr << "\n");

  try {
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, updateToSolutionVecTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    updateToSolutionVecTransStr = "Identity";

    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, updateToSolutionVecTransStr);
  }

  cf_assert(vecTransProv.isNotNull());
  _updateToSolutionVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(_updateToSolutionVecTrans.getPtr() != CFNULL);
  
  std::string solutionToLinearVecTransStr =
    VarSetTransformer::getProviderName(physModel->getConvectiveName(), _solutionVarStr, _linearVarStr);
  
  CFLog(VERBOSE, "Configuring VarSet Transformer: " << solutionToLinearVecTransStr << "\n");
  
  try {
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, solutionToLinearVecTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    solutionToLinearVecTransStr = "Identity";

    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, solutionToLinearVecTransStr);
  }
  
  cf_assert(vecTransProv.isNotNull());
  _solutionToLinearVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(_solutionToLinearVecTrans.getPtr() != CFNULL);
  
  std::string updateToReconstructVecTransStr =
    VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), _updateVarStr, _reconstructVarStr);
  
  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
	updateToReconstructVecTransStr << "\n");
  
  try {
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, updateToReconstructVecTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    updateToReconstructVecTransStr = "Identity";
    
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, updateToReconstructVecTransStr);
  }

  cf_assert(vecTransProv.isNotNull());
  _updateToReconstrVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(_updateToReconstrVecTrans.getPtr() != CFNULL);
  
  std::string reconstructToUpdateVecTransStr =
    VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), _reconstructVarStr, _updateVarStr);
  
  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
	reconstructToUpdateVecTransStr << "\n");
  
  try {
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, reconstructToUpdateVecTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    reconstructToUpdateVecTransStr = "Identity";
    
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
    vecTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, reconstructToUpdateVecTransStr);
  }

  cf_assert(vecTransProv.isNotNull());
  _reconstrToUpdateVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(_reconstrToUpdateVecTrans.getPtr() != CFNULL);
}
      
//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::configureJacobianLinearizer( Config::ConfigArgs& args )
{
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp =
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string linearizerName = physModel->getConvectiveName() + "Linear" + _linearVarStr;
  
  CFLog(VERBOSE, "CellCenterFVMData::configureJacobianLinearizer() => linearizerName = "
        << linearizerName << "\n");
  
  SafePtr<JacobianLinearizer::PROVIDER> jacobLinearProv = CFNULL;
  try {
    jacobLinearProv =
      FACTORY_GET_PROVIDER(getFactoryRegistry(), JacobianLinearizer, linearizerName);
  }
  catch (Common::NoSuchValueException& except) {
    linearizerName = "Null";
    jacobLinearProv = 
      FACTORY_GET_PROVIDER(getFactoryRegistry(), JacobianLinearizer, linearizerName);
  }
  _linearizer = jacobLinearProv->create(physModel);
  cf_assert(_linearizer.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "CellCenterFVMData::setup() => start\n");
    
  SpaceMethodData::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // set up the GeometricEntity builders
  _faceTrsGeoBuilder.setup();
  _faceCellTrsGeoBuilder.setup();
  _cellTrsGeoBuilder.setup();
  _geoWithNodesBuilder.setup();
  
  _unitNormal.resize(dim);
  
  _volumeIntegrator.setup();
  
  // sort TRSs to which BCs are not associated
  sort(_trssWithNoBC.begin(), _trssWithNoBC.end()); 
  
  CFLog(VERBOSE, "CellCenterFVMData::setup() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void CellCenterFVMData::unsetup()
{  
  SpaceMethodData::unsetup(); 
  
  _faceTrsGeoBuilder.unsetup();
  _faceCellTrsGeoBuilder.unsetup();
  _cellTrsGeoBuilder.unsetup();
  _geoWithNodesBuilder.unsetup();
  
  _volumeIntegrator.unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFuint CellCenterFVMData::getOppositeIFace(CFuint iFace, CFuint dim,
					   CFuint nbCellNodes) const
{
  if (dim == DIM_1D) {
   return (iFace == 0) ? 1 : 0;
  }
  else if (dim == DIM_2D && nbCellNodes == 4) {
    switch(iFace) {
    case 0:
      return 2;
      break;
    case 1:
      return 3;
      break;
    case 2:
      return 0;
      break;
    case 3:
      return 1;
      break;
    }
  }
  else if (dim == DIM_3D && nbCellNodes == 8) {
    switch(iFace) {
    case 0:
      return 1;
      break;
    case 1:
      return 0;
      break;
    case 2:
      return 4;
      break;
    case 3:
      return 5;
      break;
    case 4:
      return 2;
      break;
    case 5:
      return 3;
      break;
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVMData::setBCList(const vector<SelfRegistPtr<CellCenterFVMCom> >& bcList)
{
  const CFuint nbBCs = bcList.size();
  _bcMap.reserve(nbBCs);
  
  vector<FVMCC_BC*> tempBC(nbBCs);
  for (CFuint iBC = 0; iBC < nbBCs; ++iBC) {
    tempBC[iBC] = dynamic_cast<FVMCC_BC*>(bcList[iBC].getPtr());
    if(tempBC[iBC] == CFNULL){
      std::string msg = std::string("The boundary condition with the name: ") +
        bcList[iBC]->getName() +
        std::string(" is not a FVMCC boundary condition");
      throw NoSuchValueException (FromHere(),msg);
    }
  }
  
  // build a mapping iTRS -> appropriate BC
  vector<Common::SafePtr<TopologicalRegionSet> > trs =
    MeshDataStack::getActive()->getTrsList();

  const CFuint nbTRSs = trs.size();
  CFLog(VERBOSE, "FVMCC_ComputeRHS::setBCList() => nbBCs = " << nbBCs << ", nbTRSs = " << nbTRSs << "\n");
  
  for (CFuint iBC = 0; iBC < nbBCs; ++iBC) {
    // one BC can be associated to more TRSs
    // however there is only one BC associated to each TRS
    const vector<std::string>& names = tempBC[iBC]->getTrsNames();
    for (CFuint iName = 0; iName < names.size(); ++iName) {
      std::string nameTRS = names[iName];
      bool nameFound = false;
      for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
	CFLog(VERBOSE, "FVMCC_ComputeRHS::setBCList() => " << trs[iTRS]->getName() << " == " << nameTRS << "?\n");
	if (trs[iTRS]->getName() == nameTRS) {
	  CFLog(VERBOSE, "found iTRS = " << iTRS << " == " << tempBC[iBC] << "\n");
	  _bcMap.insert(iTRS, tempBC[iBC]);
	  nameFound = true;
	  break;
	}
      }
      if (!nameFound) {
	std::string msg = std::string("TRS named ") +
	  nameTRS + std::string("not found");
	throw NoSuchValueException (FromHere(),msg);
      }
    }
  }
  
  // this is fundamental when using CFMap !!!
  _bcMap.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////
 
Common::SafePtr<Common::CFMap<CFuint, FVMCC_BC*> > CellCenterFVMData::getMapBC()
{
  return &_bcMap;
} 
  
//////////////////////////////////////////////////////////////////////////////
 
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

