#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "FiniteVolume/ComputeDiffusiveFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CellCenterFVM,
			    SpaceMethod,
			    FiniteVolumeModule,
			    1>
cellCenterFVMMethodProvider("CellCenterFVM");

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
   options.addConfigOption< std::vector<std::string> >("UnSetupNames","Names of the setup commands.");
   options.addConfigOption< std::string >("ComputeTimeRHS","Compute time contribution to the right hand side.");
   options.addConfigOption< std::string >("SetNodalStatesCom","SetNodalStates Command to set the nodal states.");
   options.addConfigOption< std::string >("ComputeRHS","Compute right hand side. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("UnSetupCom","UnSetup Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
   options.addConfigOption< std::vector<std::string> >("SetupCom","Setup Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PreProcessCom","Pre-process Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("AfterMeshUpdateCom","Command to run after mesh update.");
   options.addConfigOption< std::string >("BeforeMeshUpdateCom","Command to run before mesh update.");
   options.addConfigOption< std::vector<std::string> >("BcComds","Types of the boundary conditions commands.");
   options.addConfigOption< std::vector<std::string> >("SetupNames","Names of the setup commands.");
   options.addConfigOption< std::vector<std::string> >("PreProcessNames","Names of the preprocess commands.");
   options.addConfigOption< std::vector<std::string> >("BcNames","Names for the configuration of the boundary conditions commands.");
   options.addConfigOption< bool >("OnlyInitComs","Use only init commands to initialize.");
   options.addConfigOption< std::string >("SpaceRHSForGivenCell","Command for the computation of the space discretization contibution to RHS for one cell.");
   options.addConfigOption< std::string >("TimeRHSForGivenCell" ,"Command for the computation of the space discretization contibution to RHS for one cell.");
}

//////////////////////////////////////////////////////////////////////////////

CellCenterFVM::CellCenterFVM(const std::string& name)
  : SpaceMethod(name),
    _data(new CellCenterFVMData(this)),
    _setups(),
    _unSetups(),
    _preProcess(),
    _extrapolateStates(),
    _computeSpaceRHS(),
    _computeTimeRHS(),
    _inits(0),
    _bcs(0),
    _beforeMeshUpdate(),
    _afterMeshUpdate(),
    _spaceRHSForGivenCell(),
    _timeRHSForGivenCell(),
    _isBcApplied(false)
{
  addConfigOptionsTo(this);

  cf_assert(_data.isNotNull());

  // set default value of builder for CellCenterFVM
  // to be FVMCC_MeshDataBuilder
  m_builder = "FVMCC";
  // set default global jacobian sparsity
  m_sparsity = "FVMCellCentered";

  _useOnlyInitComs = false;
  setParameter("OnlyInitComs",&_useOnlyInitComs);

  _setupStr = vector<std::string>();
  setParameter("SetupCom",&_setupStr);

  _setupNameStr = vector<std::string>();
  setParameter("SetupNames",&_setupNameStr);

  _unSetupStr = vector<std::string>();
  setParameter("UnSetupCom",&_unSetupStr);

  _unSetupNameStr = vector<std::string>();
  setParameter("UnSetupNames",&_unSetupNameStr);
  
  _unSetupStr = vector<std::string>();
  setParameter("UnSetupCom",&_unSetupStr);
  
  _unSetupNameStr = vector<std::string>();
  setParameter("UnSetupNames",&_unSetupNameStr);

   _preProcessStr = vector<std::string>();
  setParameter("PreProcessCom",&_preProcessStr);
  
  _preProcessNameStr = vector<std::string>();
  setParameter("PreProcessNames",&_preProcessNameStr);
  
  _extrapolateStatesStr = "StdSetNodalStates";
  setParameter("SetNodalStatesCom",&_extrapolateStatesStr);

  _computeSpaceRHSStr = "FVMCC";
  setParameter("ComputeRHS",&_computeSpaceRHSStr);

  _computeTimeRHSStr = "Null";
  setParameter("ComputeTimeRHS",&_computeTimeRHSStr);

  _initTypeStr.clear();
  setParameter("InitComds",&_initTypeStr);

  _initNameStr.clear();
  setParameter("InitNames",&_initNameStr);

  _bcTypeStr.clear();
  setParameter("BcComds",&_bcTypeStr);

  _bcNameStr.clear();
  setParameter("BcNames",&_bcNameStr);

  // Moving Mesh
  _beforeMeshUpdateStr = "Null";
  setParameter("BeforeMeshUpdateCom",&_beforeMeshUpdateStr);

  _afterMeshUpdateStr = "Null";
  setParameter("AfterMeshUpdateCom",&_afterMeshUpdateStr);

  // default values for LU-SGS-related commands
  setParameter( "SpaceRHSForGivenCell", &_spaceRHSForGivenCellStr);
  _spaceRHSForGivenCellStr = "Null";

  setParameter( "TimeRHSForGivenCell", &_timeRHSForGivenCellStr);
  _timeRHSForGivenCellStr = "Null";
}

//////////////////////////////////////////////////////////////////////////////

CellCenterFVM::~CellCenterFVM()
{
  CFAUTOTRACE;

  clearSetupComs();
  clearUnSetupComs();
  clearPreProcessComs();
  clearInitComs();
  clearBCComs();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> CellCenterFVM::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<SpaceMethodData> CellCenterFVM::getSpaceMethodData()
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::setCollaborator(MultiMethodHandle<LinearSystemSolver> lss)
{
  _data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  _data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::clearInitComs()
{
  CFAUTOTRACE;

   vector<SelfRegistPtr<CellCenterFVMCom> >().swap(_inits);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::clearBCComs()
{
  CFAUTOTRACE;

  vector<SelfRegistPtr<CellCenterFVMCom> >().swap(_bcs);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::clearSetupComs()
{
  CFAUTOTRACE;

  vector<SelfRegistPtr<CellCenterFVMCom> >().swap(_setups);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::clearUnSetupComs()
{
  CFAUTOTRACE;

  vector<SelfRegistPtr<CellCenterFVMCom> >().swap(_unSetups);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::clearPreProcessComs()
{
  CFAUTOTRACE;

  vector<SelfRegistPtr<CellCenterFVMCom> >().swap(_preProcess);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "CellCenterFVM::configure()\n");
  SpaceMethod::configure(args);
  
  CFLog(VERBOSE, "CellCenterFVM::configureNested()\n");
  
  _data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _data.getPtr(), args );
  
  // add here configures to the CellCenterFVM
  clearSetupComs();
  clearUnSetupComs();
  clearPreProcessComs();
  
  cf_assert(_setupStr.size() == _setupNameStr.size());
  if(_setupStr.size() == 0){
    _setupStr.resize(1);
    _setupNameStr.resize(1);
    _setupStr[0] = "StdSetup";
    _setupNameStr[0] = "StdSetup1";
  }
  
  _setups.resize(_setupStr.size());

  for(CFuint i = 0; i < _setups.size(); ++i) {

    CFLog(INFO, "SETUP type = " << _setupStr[i] << "\n");
    CFLog(INFO, "SETUP name = " << _setupNameStr[i] << "\n");

    configureCommand<CellCenterFVMCom,
      CellCenterFVMData,
      CellCenterFVMComProvider>
      (args, _setups[i], _setupStr[i],_setupNameStr[i], _data);

    cf_assert(_setups[i].isNotNull());
  }

  cf_assert(_unSetupStr.size() == _unSetupNameStr.size());
  if(_setupStr.size() == 0){
    _unSetupStr.resize(1);
    _unSetupNameStr.resize(1);
    _unSetupStr[0] = "StdUnSetup";
    _unSetupNameStr[0] = "StdUnSetup1";
  }

  _unSetups.resize(_unSetupStr.size());

  for(CFuint i = 0; i < _unSetups.size(); ++i) {

    CFLog(INFO, "SETUP type = " << _unSetupStr[i] << "\n");
    CFLog(INFO, "SETUP name = " << _unSetupNameStr[i] << "\n");

    configureCommand<CellCenterFVMCom,
      CellCenterFVMData,
      CellCenterFVMComProvider>
      (args, _unSetups[i], _unSetupStr[i],_unSetupNameStr[i], _data);

    cf_assert(_unSetups[i].isNotNull());
  }

  cf_assert(_preProcessStr.size() == _preProcessNameStr.size());
  if(_preProcessStr.size() == 0){
    _preProcessStr.resize(1);
    _preProcessNameStr.resize(1);
    _preProcessStr[0] = "Null";
    _preProcessNameStr[0] = "PreProcess1";
  }
  
  _preProcess.resize(_preProcessStr.size());

  for(CFuint i = 0; i < _preProcess.size(); ++i) {

    CFLog(INFO, "PREPROCESS type = " << _preProcessStr[i] << "\n");
    CFLog(INFO, "PREPROCESS name = " << _preProcessNameStr[i] << "\n");

    configureCommand<CellCenterFVMCom,
      CellCenterFVMData,
      CellCenterFVMComProvider>
      (args, _preProcess[i], _preProcessStr[i],_preProcessNameStr[i], _data);

    cf_assert(_preProcess[i].isNotNull());
  }
  
  
  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _extrapolateStates,
                                                              _extrapolateStatesStr,
                                                              _data);

  CFLog(INFO,"CellCenterFVM: Using ComputeRHS: " << _computeSpaceRHSStr << "\n");
  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _computeSpaceRHS,
							       _computeSpaceRHSStr,
							       _data);

  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _computeTimeRHS,
							       _computeTimeRHSStr,
							       _data);

  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _beforeMeshUpdate,
							       _beforeMeshUpdateStr,
							       _data);

  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _afterMeshUpdate,
							       _afterMeshUpdateStr,
							       _data);

  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _spaceRHSForGivenCell,
                                                               _spaceRHSForGivenCellStr,
                                                               _data);

  configureCommand<CellCenterFVMData,CellCenterFVMComProvider>(args, _timeRHSForGivenCell,
                                                               _timeRHSForGivenCellStr,
                                                               _data);

  clearInitComs();
  clearBCComs();

  cf_assert(_initTypeStr.size() == _initNameStr.size());

  _inits.resize(_initTypeStr.size());
  for(CFuint i = 0; i < _inits.size(); ++i) {

    CFLog(INFO, "INIT type = " << _initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << _initNameStr[i] << "\n");

    configureCommand<CellCenterFVMCom,
      CellCenterFVMData,
      CellCenterFVMComProvider>(args,
                                _inits[i],
                                _initTypeStr[i],
                                _initNameStr[i],
                                _data);
    cf_assert(_inits[i].isNotNull());
  }

  cf_assert(_bcTypeStr.size() == _bcNameStr.size());

  _bcs.resize(_bcTypeStr.size());
  if (_bcs.size() == 0) {
    _bcTypeStr.push_back("NullBC");
    _bcNameStr.push_back("DummyBC");
    _bcs.resize(_bcTypeStr.size());
  }

  for(CFuint i = 0; i < _bcs.size(); ++i) {

    CFLog(INFO, "BC type = " << _bcTypeStr[i] << "\n");
    CFLog(INFO, "BC name = " << _bcNameStr[i] << "\n");

    configureCommand<CellCenterFVMCom,
      CellCenterFVMData,
      CellCenterFVMComProvider>(args,
                                _bcs[i],
                                _bcTypeStr[i],
                                _bcNameStr[i],
                                _data);
    cf_assert(_bcs[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;
  
  _data->setIsInitializationPhase(true);
    
  // always initialize if it is not a restart
  if (!isRestart) {
    for(CFuint i = 0; i < _inits.size(); ++i) {
      cf_assert(_inits[i].isNotNull());
      CFLog(VERBOSE, "Initializing " << _inits[i]->getName() << "START\n");
      _inits[i]->execute();
      CFLog(VERBOSE, "Initializing " << _inits[i]->getName() << "END\n");
    }
    
    if (!_useOnlyInitComs) {
      applyBCImpl();
    }
  }
  else {
    if (!_useOnlyInitComs) {
      applyBCImpl();
    }
    else {
      // if it is a restart, don't initialize the internal State's
      for(CFuint i = 0; i < _inits.size(); ++i) {
	cf_assert(_inits[i].isNotNull());
	
	// check the names of the TRS on which to apply the init Command
	const vector<std::string>& names = _inits[i]->getTrsNames();
	vector<std::string>::const_iterator it;
	bool toInitialize = true;
	for (it = names.begin(); it != names.end(); ++it) {
	  if (*it == "InnerFaces" || *it == "InnerCells") {
	    toInitialize = false;
	    break;
	  }
	}
	
	if (toInitialize) {
	  _inits[i]->execute();
	}
      }
    }
  }
  _data->setIsInitializationPhase(false);
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
  
  if (!_data->doOnlyPreprocessSolution()) {
    // preprocess solution
    preProcessSolution();
  }
    
  // set the residual factor for which residual and jacobian have to
  // be multiplied
  _data->setResFactor(factor);
  
  // apply the BC
  applyBC();
  
  // BC should actually be applied after the computeResidual
  // and after the update of the states !!!
  _computeSpaceRHS->execute();
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  // set the residual factor for which residual and jacobian have to
  // be multiplied
  _data->setResFactor(factor);

  if (!_computeTimeRHS->isNull()) {
    _computeTimeRHS->execute();
  }
  checkMatrixFrozen();
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::applyBCImpl()
{
  CFAUTOTRACE;

  cf_assert(isSetup());
  cf_assert(isConfigured());
  
  for(CFuint i = 0; i < _bcs.size(); ++i) {
    cf_assert(_bcs[i].isNotNull());
    CFLog(VERBOSE, "Applying BC " << _bcs[i]->getName() << " => START\n");
    _bcs[i]->execute();
    CFLog(VERBOSE, "Applying BC " << _bcs[i]->getName() << " => END\n");
  }
  
  // reset to false the preProcessBCFlag
  _data->preProcessBCFlag() = false;
  _isBcApplied = true;
}
      
//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::prepareComputationImpl()
{
  CFAUTOTRACE;
  // currently doing nothing
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::postProcessSolutionImpl()
{
  CFAUTOTRACE;
  // currently doing nothing
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::preProcessSolutionImpl()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "CellCenterFVM::preProcessSolutionImpl() => START\n");
  
  _data->setIsPreProcessedSolution(false);
  
  // pre-process commands
  for(CFuint i=0; i < _preProcess.size();i++){
    cf_assert(_preProcess[i].isNotNull());
    CFLog(VERBOSE, "CellCenterFVM::computeSpaceResidualImpl() => pre-processing start"
	  << _preProcess[i]->getName() << " \n");
    _preProcess[i]->execute();
    CFLog(VERBOSE, "CellCenterFVM::computeSpaceResidualImpl() => pre-processing end"
	  << _preProcess[i]->getName() << " \n");
  }

  _data->setIsPreProcessedSolution(true);
  
  CFLog(VERBOSE, "CellCenterFVM::preProcessSolutionImpl() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::setMethodImpl()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "CellCenterFVM::setMethodImpl() => before SpaceMethod::setMethodImpl()\n");
  
  SpaceMethod::setMethodImpl();
  
  CFLog(VERBOSE, "CellCenterFVM::setMethodImpl() => after SpaceMethod::setMethodImpl()\n");
 
  //setup has to be processed first to create InwardNormals
  // needed by the BC commands to create the nodal normals
  // in Method::setupCommandsAndStrategies()
  
  // store the mapping index TRS -> bc commands into the data
  _data->setBCList(_bcs);
  
  for(CFuint i=0; i < _setups.size();i++){
    cf_assert(_setups[i].isNotNull());
    CFLog(VERBOSE, "CellCenterFVM::setMethodImpl() => start setting up " << _setups[i]->getName() << " \n");
    _setups[i]->execute();
    CFLog(VERBOSE, "CellCenterFVM::setMethodImpl() => end setting up " << _setups[i]->getName() << " \n");
  }
  
  CFLog(VERBOSE, "CellCenterFVM::setMethodImpl() => before setupCommandsAndStrategies()\n");
  setupCommandsAndStrategies();
  CFLog(VERBOSE, "CellCenterFVM::setMethodImpl() => after setupCommandsAndStrategies()\n");
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::unsetMethodImpl()
{
  CFAUTOTRACE;

  unsetupCommandsAndStrategies();

  for(CFuint i=0; i < _unSetups.size();++i){
    cf_assert(_unSetups[i].isNotNull());
    _unSetups[i]->execute();
  }

  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t CellCenterFVM::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;

  _beforeMeshUpdate->execute();

  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t CellCenterFVM::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;

  _afterMeshUpdate->execute();

  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;

  // ghost states have to be updated before extrapolating to nodes for output
  if (!_isBcApplied) {applyBCImpl();}
  // applyBCImpl();
  _extrapolateStates->execute();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > CellCenterFVM::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result;

  // add strategies here
  result.push_back(_data->getPolyReconstructor().d_castTo<NumericalStrategy>());
  result.push_back(_data->getLimiter().d_castTo<NumericalStrategy>());
  result.push_back(_data->getNodalStatesExtrapolator().d_castTo<NumericalStrategy>());
  result.push_back(_data->getFluxSplitter().d_castTo<NumericalStrategy>());
  result.push_back(_data->getGeoDataComputer().d_castTo<NumericalStrategy>());
  result.push_back(_data->getDerivativeComputer().d_castTo<NumericalStrategy>());
  result.push_back(_data->getDiffusiveFluxComputer().d_castTo<NumericalStrategy>());
  
  SafePtr<vector<SelfRegistPtr<ComputeSourceTerm<CellCenterFVMData> > > > sourceTerms =
    _data->getSourceTermComputer();
  
  for(CFuint i=0; i<sourceTerms->size();++i){
    SafePtr<ComputeSourceTerm<CellCenterFVMData> > sourceTerm = ((*sourceTerms)[i]).getPtr();
    result.push_back(sourceTerm.d_castTo<NumericalStrategy>());
  }
  
  SafePtr<vector<SelfRegistPtr<EquationFilter<CellCenterFVMData> > > > equationFilters =
    _data->getEquationFilters();
  
  for(CFuint i=0; i<equationFilters->size();++i){
    SafePtr<EquationFilter<CellCenterFVMData> > eqFilter = ((*equationFilters)[i]).getPtr();
    result.push_back(eqFilter.d_castTo<NumericalStrategy>());
  }
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::checkMatrixFrozen() const
{
  if(!_data->isSysMatrixFrozen() && _data->isSysMatrixFreezedEveryIteration()) {
    CFLog(INFO,"Freezing System Matrix" << "\n");
    _data->setSysMatrixFrozen(true);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::computeSpaceRhsForStatesSetImpl(CFreal factor)
{
  // set the residual factor in the MethodData
  _data->setResFactor(factor);

  cf_assert(_spaceRHSForGivenCell.isNotNull());
  _spaceRHSForGivenCell->execute();
}

//////////////////////////////////////////////////////////////////////////////

void CellCenterFVM::computeTimeRhsForStatesSetImpl(CFreal factor)
{ 
  // set the residual factor in the MethodData
  _data->setResFactor(factor);

  cf_assert(_timeRHSForGivenCell.isNotNull());
  _timeRHSForGivenCell->execute();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
