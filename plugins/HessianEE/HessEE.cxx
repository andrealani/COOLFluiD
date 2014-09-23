#include "HessEE.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "HessianEE/HessianEE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  using namespace Framework;

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<HessEE,
               ErrorEstimatorMethod,
               HessianEEModule,
               1>
hessianErrorEstimatorProvider("HessianEE");

//////////////////////////////////////////////////////////////////////////////

void HessEE::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UpdateCom","Command to update global metric field.");
   options.addConfigOption< std::string >("FunctionCom","Command to compute the adapted function.");
   options.addConfigOption< std::string >("MetricCom","Command to compute the Metric.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("HessianCom","Command to compute the Hessian.");
   options.addConfigOption< std::string >("SmootherCom","Command to smooth either the Hessian or the Metric.");
}

//////////////////////////////////////////////////////////////////////////////

HessEE::HessEE(const std::string& name)
  : ErrorEstimatorMethod(name)
{
   addConfigOptionsTo(this);

  _data.reset( new HessEEData(this) );

  _setupStr = "StdSetup";
   setParameter("SetupCom",&_setupStr);

  _unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&_unSetupStr);

  _computeFunctionStr = "FunctionCalc";
   setParameter("FunctionCom",&_computeHessianStr);

  _computeHessianStr = "IntegralHessCalc";
   setParameter("HessianCom",&_computeHessianStr);

  _computeMetricStr = "StdMetricCalc";
   setParameter("MetricCom",&_computeMetricStr);

  _smootherStr = "StdSmoothCom";
   setParameter("SmootherCom",&_smootherStr);

  _updaterStr = "StdUpdateCom";
   setParameter("UpdateCom",&_updaterStr);
}

//////////////////////////////////////////////////////////////////////////////

HessEE::~HessEE()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> HessEE::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void HessEE::configure ( Config::ConfigArgs& args )
{
  ErrorEstimatorMethod::configure(args);
  configureNested ( _data.getPtr(), args );

  // add configures to the HessEECom's

  configureCommand<HessEEData,HessEEComProvider>( args, _setup,_setupStr,_data);

  configureCommand<HessEEData,HessEEComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<HessEEData,HessEEComProvider>( args, _computeFunction,_computeFunctionStr,_data);

  configureCommand<HessEEData,HessEEComProvider>( args, _computeHessian,_computeHessianStr,_data);

  configureCommand<HessEEData,HessEEComProvider>( args, _computeMetric,_computeMetricStr,_data);

  configureCommand<HessEEData,HessEEComProvider>( args, _smoother,_smootherStr,_data);

  configureCommand<HessEEData,HessEEComProvider>( args, _updater,_updaterStr,_data);
}

//////////////////////////////////////////////////////////////////////////////

void HessEE::setMethodImpl()
{
  ErrorEstimatorMethod::setMethodImpl();

  _data->getGeoWithNodesBuilder()->setup();
  setupCommandsAndStrategies();
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void HessEE::unsetMethodImpl()
{
  _unSetup->execute();
  unsetupCommandsAndStrategies();

  ErrorEstimatorMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void HessEE::estimateImpl()
{

  CFAUTOTRACE;

  if(!(SubSystemStatusStack::getActive()->getNbIter() % _estimateRate))
  {
    // set the nodal states at first
    getMethodData()->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();

    _computeFunction->execute();

    _computeHessian->execute();

    _data->rSmthDataName() = "hessian";
    _data->rSmthNIter() = _data->getHessSmoothNb();
    _data->rSmthWght() = _data->getHessSmoothWght();
    _smoother->execute();

    _computeMetric->execute();

    _data->rSmthDataName() = "metric";
    _data->rSmthNIter() = _data->getMetricSmoothNb();
    _data->rSmthWght() = _data->getMetricSmoothWght();
    _smoother->execute();

    _updater->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

