#include "HessianEE/HessianEE.hh"

#include "HessEEData.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<HessEEData>, HessEEData, HessianEEModule> nullHessEEComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void HessEEData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("InvertMetricSmooth","Invert metric for smoothing?");
   options.addConfigOption< CFuint >("NbMetricSmooth","Number of iterations of metric smoothing");
   options.addConfigOption< bool >("InvertHessSmooth","Invert hessian for smoothing?");
   options.addConfigOption< CFreal >("MaxMetricAR","Maximum Aspect Ratio for limiting the metric");
   options.addConfigOption< CFreal >("MinMetricLimit","Minimum value for limiting the metric");
   options.addConfigOption< CFreal >("Constant","Magic constant");
   options.addConfigOption< CFreal >("WghtMetricSmooth","Weight used in metric smoothing");
   options.addConfigOption< CFreal >("WghtHessSmooth","Weight used in hessian smoothing");
   options.addConfigOption< bool >("SmoothHessian","Will we smooth the Hessian?");
   options.addConfigOption< CFreal >("MaxMetricLimit","Maximum value for limiting the metric");
   options.addConfigOption< CFuint >("NbHessSmooth","Number of iterations of hessian smoothing");
}

//////////////////////////////////////////////////////////////////////////////

HessEEData::HessEEData(Common::SafePtr<Framework::Method> owner)
 : ErrorEstimatorData(owner),
   _geoWithNodesBuilder()
{
   addConfigOptionsTo(this);
  _isHessianSmooth = true;
   setParameter("SmoothHessian",&_isHessianSmooth);

  _metricLimitMax = 1.0;
   setParameter("MaxMetricLimit",&_metricLimitMax);

  _metricLimitMin = 0.0;
   setParameter("MinMetricLimit",&_metricLimitMin);

  _metricLimitAR = 10.0;
   setParameter("MaxMetricAR",&_metricLimitAR);

  _hessSmthNIter = 0;
   setParameter("NbHessSmooth",&_hessSmthNIter);

  _hessSmthWght = 0.05;
   setParameter("WghtHessSmooth",&_hessSmthWght);

  _bHessSmthInvert = false;
   setParameter("InvertHessSmooth",&_bHessSmthInvert);

  _metricSmthNIter = 0;
   setParameter("NbMetricSmooth",&_metricSmthNIter);

  _metricSmthWght = 0.05;
   setParameter("WghtMetricSmooth",&_metricSmthWght);

  _bMetricSmthInvert = false;
   setParameter("InvertMetricSmooth",&_bMetricSmthInvert);

  _magConst = 2.5;
   setParameter("Constant",&_magConst);
}

//////////////////////////////////////////////////////////////////////////////

HessEEData::~HessEEData()
{
}

//////////////////////////////////////////////////////////////////////////////

void HessEEData::configure ( Config::ConfigArgs& args )
{
  ErrorEstimatorData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

