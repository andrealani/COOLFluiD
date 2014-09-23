#ifndef COOLFluiD_Numerics_HessianEE_HessEEData_hh
#define COOLFluiD_Numerics_HessianEE_HessEEData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/ErrorEstimatorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * HessianEECom 's that compose the HessianEE.
   *
   * @see HessianEECom
   *
   * @author Tiago Quintino
   * @author Jurek Majewski
   */
class HessEEData : public Framework::ErrorEstimatorData {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  HessEEData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~HessEEData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "HessEE";
  }

  /**
   * @return if hessian should be smoothed
   */
  bool isHessianSmooth()
  {
    return _isHessianSmooth;
  }

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> >
  getGeoWithNodesBuilder()
  {
    return &_geoWithNodesBuilder;
  }

  /**
   * @return maximum value to which the metric should be limited
   */
  const CFreal& getMetricMaxLimit() const	{ return _metricLimitMax; }

  /**
   * @return minimum value to which the metric should be limited
   */
  const CFreal& getMetricMinLimit() const	{ return _metricLimitMin; }

  /**
   * @return maximum aspect ratio of a cell to which the metric should be limited
   */
  const CFreal& getMetricMaxAR() const		{ return _metricLimitAR; }

  /**
   * @return number of metric smoothing iterations
   */
  const CFuint& getMetricSmoothNb() const 	{ return _metricSmthNIter; }

  /**
   * @return weight used for metric smoothing iterations
   */
  const CFreal& getMetricSmoothWght() const	{ return _metricSmthWght; }

  /**
   * @return if the metric should be inverted before smoothing
   */
  const bool& getMetricInvert() const 	{ return _bMetricSmthInvert; }

  /**
   * @return number of hessian smoothing iterations
   */
  const CFuint& getHessSmoothNb() const 	{ return _hessSmthNIter; }

  /**
   * @return weight used for hessian smoothing iterations
   */
  const CFreal& getHessSmoothWght() const	{ return _hessSmthWght; }

  /**
   * @return if the hessian should be inverted before smoothing
   */
  const bool& getHessInvert() const 		{ return _bHessSmthInvert; }

  /**
   * @return magic constant
   */
  const CFreal& getConstant() const	{ return _magConst; }

  std::string& rSmthDataName()	    { return _smthDataName; }

  CFuint&   rSmthNIter()	    { return _smthNIter; }

  CFreal&   rSmthWght()	        { return _smthWght; }

private:

  // builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _geoWithNodesBuilder;

  /// a flag to check if we smooth the Hessian after the computation
  bool _isHessianSmooth;

  CFreal	_metricLimitMax;	// metric limit. - max spacing
  CFreal	_metricLimitMin;	// metric limit. - min spacing
  CFreal	_metricLimitAR;		// metric limit. - max aspect ratio

  CFuint	_metricSmthNIter;	// metric smoothing: no of iterations
  CFreal	_metricSmthWght;	// metric smoothing: weight
  bool	_bMetricSmthInvert;	// metric smoothing: inverted or not

  CFuint	_hessSmthNIter;		// hessian smoothing: no of iterations
  CFreal	_hessSmthWght;		// hessian smoothing: weight
  bool	_bHessSmthInvert;	// hessian smoothing: inverted or not

  CFreal	_magConst;			// magic constant

  std::string      _smthDataName;
  CFuint	_smthNIter;
  CFreal	_smthWght;

}; // end of class HessEEData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for HessianEE
typedef Framework::MethodCommand<HessEEData> HessEECom;

/// Definition of a command provider for HessianEE
typedef Framework::MethodCommand<HessEEData>::PROVIDER HessEEComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_HessEEData_hh
