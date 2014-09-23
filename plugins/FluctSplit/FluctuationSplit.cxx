#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"

#include "FluctSplit/FluctuationSplit.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/ComputeJacobStrategy.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"
#include "FluctSplit/ComputeSourceTermFSM.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FluctuationSplit,
               SpaceMethod,
               FluctSplitModule,
               1>
fluctuationSplitSpaceMethodProvider("FluctuationSplit");

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("UnSetupNames","Names of the setup commands.");
  options.addConfigOption< std::string >("ComputeTimeRHS","Compute time contribution to the right hand side.");
  options.addConfigOption< std::string >("SetNodalStatesCom","SetNodalStates Command to set the nodal states.");
  options.addConfigOption< std::string >("ComputeRHS","Compute right hand side. This command seldomly needs overriding.");
  options.addConfigOption< std::vector<std::string> >("UnSetupCom","UnSetup Command to run. This command seldomly needs overriding.");
  options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("SetupCom","Setup Command to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("AfterMeshUpdateCom","Command to run after mesh update.");
  options.addConfigOption< std::string >("BeforeMeshUpdateCom","Command to run before mesh update.");
  options.addConfigOption< std::vector<std::string> >("BcComds","Types of the boundary conditions commands.");
  options.addConfigOption< std::vector<std::string> >("SetupNames","Names of the setup commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names for the configuration of the boundary conditions commands.");
}

//////////////////////////////////////////////////////////////////////////////

FluctuationSplit::FluctuationSplit(const std::string& name)
  : SpaceMethod(name),
    _setups(),
    _unSetups(),
    _extrapolateStates(),
    _computeSpaceRHS(),
    _computeTimeRHS(),
    _inits(0),
    _bcs(0)
{
  addConfigOptionsTo(this);
  _data.reset(new FluctuationSplitData(this));
  cf_assert(_data.isNotNull());

  // set default value of builder for FlucutationSplit
  // to be RDS_MeshDataBuilder
  m_builder = "RDS";
  // set default global jacobian sparsity
  m_sparsity = "CellVertex";

  _setupStr = vector<std::string>();
  setParameter("SetupCom",&_setupStr);

  _setupNameStr = vector<std::string>();
  setParameter("SetupNames",&_setupNameStr);

  _unSetupStr = vector<std::string>();
  setParameter("UnSetupCom",&_unSetupStr);

  _unSetupNameStr = vector<std::string>();
  setParameter("UnSetupNames",&_unSetupNameStr);

  _extrapolateStatesStr = "Null";
  setParameter("SetNodalStatesCom",&_extrapolateStatesStr);

  _computeSpaceRHSStr = "StdComputeRHS";
  setParameter("ComputeRHS",&_computeSpaceRHSStr);

  _computeTimeRHSStr = "Null";
  setParameter("ComputeTimeRHS",&_computeTimeRHSStr);

  _initTypeStr = vector<std::string>();
  setParameter("InitComds",&_initTypeStr);

  _initNameStr = vector<std::string>();
  setParameter("InitNames",&_initNameStr);

  _bcTypeStr = vector<std::string>();
  setParameter("BcComds",&_bcTypeStr);

  _bcNameStr = vector<std::string>();
  setParameter("BcNames",&_bcNameStr);

  // Moving Mesh
  _beforeMeshUpdateStr = "Null";
  setParameter("BeforeMeshUpdateCom",&_beforeMeshUpdateStr);

  _afterMeshUpdateStr = "Null";
  setParameter("AfterMeshUpdateCom",&_afterMeshUpdateStr);

}

//////////////////////////////////////////////////////////////////////////////

FluctuationSplit::~FluctuationSplit()
{
  clearInitComs();
  clearBCComs();
  clearSetupComs();
  clearUnSetupComs();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::setCollaborator(MultiMethodHandle<LinearSystemSolver> lss)
{
  CFAUTOTRACE;
  _data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  CFAUTOTRACE;
  _data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::clearInitComs()
{
  vector<SelfRegistPtr<FluctuationSplitCom> >().swap(_inits);
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::clearBCComs()
{
  vector<SelfRegistPtr<FluctuationSplitCom> >().swap(_bcs);
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::clearSetupComs()
{
  vector<SelfRegistPtr<FluctuationSplitCom> >().swap(_setups);
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::clearUnSetupComs()
{
  vector<SelfRegistPtr<FluctuationSplitCom> >().swap(_unSetups);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> FluctuationSplit::getMethodData () const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<SpaceMethodData> FluctuationSplit::getSpaceMethodData()
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::configureSetupComs( Config::ConfigArgs& args )
{
  clearSetupComs();

  // create the standard setup command if noe is chosen
  if( _setupStr.empty() )
  {
    _setupStr.resize(1);
    _setupStr[0] = "StdSetup";
  }

  // assign names is names not given
  if ( _setupNameStr.empty() )
  {
    for ( CFuint i = 0 ; i < _setupStr.size() ; ++i )
    {
      _setupNameStr.push_back( _setupStr[i] + Common::StringOps::to_str(i) );
    }
  }

  // signal error if sizes dont match
  if (_setupStr.size() != _setupNameStr.size())
    throw BadValueException (FromHere(), "Number of names for setup commands do not match the number of types" );

  _setups.resize(_setupStr.size());

  // do the actual setup
  for(CFuint i = 0; i < _setups.size(); ++i)
  {

    CFLog(VERBOSE, "SETUP type = " << _setupStr[i] << "\n");
    CFLog(VERBOSE, "SETUP name = " << _setupNameStr[i] << "\n");

    configureCommand<FluctuationSplitCom,
      FluctuationSplitData,
      FluctuationSplitComProvider>
      (args, _setups[i], _setupStr[i],_setupNameStr[i], _data);

    cf_assert(_setups[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::configureUnsetupComs(Config::ConfigArgs& args)
{
  clearUnSetupComs();

  // create the standard ununSetup command if noe is chosen
  if( _unSetupStr.empty() )
  {
    _unSetupStr.resize(1);
    _unSetupStr[0] = "StdUnSetup";
  }

  // assign names is names not given
  if ( _unSetupNameStr.empty() )
  {
    for ( CFuint i = 0 ; i < _unSetupStr.size() ; ++i )
    {
      _unSetupNameStr.push_back( _unSetupStr[i] + Common::StringOps::to_str(i) );
    }
  }

  // signal error if sizes dont match
  if (_unSetupStr.size() != _unSetupNameStr.size())
    throw BadValueException (FromHere(), "Number of names for Unsetup commands do not match the number of types" );

  _unSetups.resize(_unSetupStr.size());

  // do the actual ununSetup
  for(CFuint i = 0; i < _unSetups.size(); ++i)
  {

    CFLog(VERBOSE, "Unsetup type = " << _unSetupStr[i] << "\n");
    CFLog(VERBOSE, "Unsetup name = " << _unSetupNameStr[i] << "\n");

    configureCommand<FluctuationSplitCom,
      FluctuationSplitData,
      FluctuationSplitComProvider>
      (args, _unSetups[i], _unSetupStr[i],_unSetupNameStr[i], _data);

    cf_assert(_unSetups[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  SpaceMethod::configure(args);
  
  configureNested ( _data.getPtr(), args );
  
  configureSetupComs(args);
  configureUnsetupComs(args);

  configureCommand<FluctuationSplitData, FluctuationSplitComProvider>  (args, _extrapolateStates, _extrapolateStatesStr, _data);
  configureCommand<FluctuationSplitData, FluctuationSplitComProvider>  (args, _computeSpaceRHS,   _computeSpaceRHSStr,   _data);
  configureCommand<FluctuationSplitData, FluctuationSplitComProvider>  (args, _computeTimeRHS,    _computeTimeRHSStr,    _data);
  configureCommand<FluctuationSplitData, FluctuationSplitComProvider>  (args, _beforeMeshUpdate,  _beforeMeshUpdateStr,  _data);
  configureCommand<FluctuationSplitData, FluctuationSplitComProvider>  (args, _afterMeshUpdate,   _afterMeshUpdateStr,   _data);

  clearInitComs();

  if (_initTypeStr.size() != _initNameStr.size())
    throw BadValueException (FromHere(), "Number of names for Init commands do not match the number of types" );

  _inits.resize(_initTypeStr.size());
  for(CFuint i = 0; i < _inits.size(); ++i) {

    CFLog(INFO, "INIT type = " << _initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << _initNameStr[i] << "\n");

    configureCommand<FluctuationSplitCom,
      FluctuationSplitData,
      FluctuationSplitComProvider>
      (args, _inits[i], _initTypeStr[i],_initNameStr[i], _data);

    cf_assert(_inits[i].isNotNull());
  }

  clearBCComs();
  if (_bcTypeStr.size() != _bcNameStr.size())
    throw BadValueException (FromHere(), "Number of names for BC commands do not match the number of types" );

  _bcs.resize(_bcTypeStr.size());
  for(CFuint i = 0; i < _bcs.size(); ++i) {

    CFLog(INFO, "BC type = " << _bcTypeStr[i] << "\n");
    CFLog(INFO, "BC name = " << _bcNameStr[i] << "\n");

    configureCommand<FluctuationSplitCom,
      FluctuationSplitData,
      FluctuationSplitComProvider>
      (args, _bcs[i],_bcTypeStr[i], _bcNameStr[i], _data);

    cf_assert(_bcs[i].isNotNull());
  } 
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;

  if (!isRestart) {
    _data->setInitializationPhase(true);
    for(CFuint i = 0; i < _inits.size(); ++i) {
      cf_assert(_inits[i].isNotNull());
      _inits[i]->execute();
    }
    _data->setInitializationPhase(false);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;

  _extrapolateStates->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  // set the residual factor for which residual and jacobian have to be multiplied
  _data->setResFactor(factor);

  cf_assert(_computeSpaceRHS.isNotNull());
  _computeSpaceRHS->execute();
  // apply the boundary conditions
  applyBC();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  // set the residual factor for which residual and jacobian have to be multiplied
  _data->setResFactor(factor);

  cf_assert(_computeTimeRHS.isNotNull());

  // if there time discretization is defined
  if (!_computeTimeRHS->isNull())
  {
    _computeTimeRHS->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::applyBCImpl()
{
  CFAUTOTRACE;

  for(CFuint i = 0; i < _bcs.size(); ++i)
  {
    cf_assert(_bcs[i].isNotNull());
    _bcs[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::prepareComputationImpl()
{
  CFAUTOTRACE;

  getData()->getFluctSplitStrategy()->prepare();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::setMethodImpl()
{
  CFAUTOTRACE;

  SpaceMethod::setMethodImpl();

  // assign all the Inner TRS's to the commands that compute RHS and Jacob
  std::vector<std::string> tags;
  tags.push_back ( "inner" );
  tags.push_back ( "cell" );
  std::vector< Common::SafePtr<TopologicalRegionSet> > inner_cell_trs = MeshDataStack::getActive()->getFilteredTrsList(tags);

  _computeSpaceRHS->setTrsList( inner_cell_trs );
  _computeTimeRHS->setTrsList( inner_cell_trs );

  // setup commands have to be processed first to create data
  // that is used by the strategies (example the InwardNormals)
  for(CFuint i=0; i < _setups.size(); ++i)
  {
    cf_assert(_setups[i].isNotNull());
    _setups[i]->execute();
  }

  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplit::unsetMethodImpl()
{
  CFAUTOTRACE;

  // reverse order that of the setup
  unsetupCommandsAndStrategies();

  for(CFuint i=0; i < _unSetups.size(); ++i){
    cf_assert(_unSetups[i].isNotNull());
    _unSetups[i]->execute();
  }

  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t FluctuationSplit::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;

  cf_assert(_beforeMeshUpdate.isNotNull());
  _beforeMeshUpdate->execute();

  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t FluctuationSplit::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;

  cf_assert(_afterMeshUpdate.isNotNull());
  _afterMeshUpdate->execute();

  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::NumericalStrategy> > FluctuationSplit::getStrategyList() const
{
  vector<Common::SafePtr<Framework::NumericalStrategy> > result;

  // add strategies here
  result.push_back(_data->getFluctSplitStrategy().d_castTo<NumericalStrategy>());
  result.push_back(_data->getArtificialDiffusionStrategy().d_castTo<NumericalStrategy>());
  result.push_back(_data->getJacobStrategy().d_castTo<NumericalStrategy>());
  result.push_back(_data->getDiffusiveTermComputer().d_castTo<NumericalStrategy>());

  const CFuint stsize = _data->getSourceTermComputer()->size();
  for (CFuint i = 0; i < stsize; ++i) {
    SafePtr<ComputeSourceTermFSM> st =
      (*_data->getSourceTermComputer())[i].getPtr();
    result.push_back(st.d_castTo<NumericalStrategy>());
  }

  result.push_back(_data->getSysSplitter().d_castTo<NumericalStrategy>());
  result.push_back(_data->getScalarSplitter().d_castTo<NumericalStrategy>());

  for (CFuint i = 0; i < stsize; ++i) {
    SafePtr<SourceTermSplitter> st =
      (*_data->getSourceTermSplitter())[i].getPtr();
    result.push_back(st.d_castTo<NumericalStrategy>());
  }

  result.push_back(_data->getJacobianFixComputer().d_castTo<NumericalStrategy>());

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
