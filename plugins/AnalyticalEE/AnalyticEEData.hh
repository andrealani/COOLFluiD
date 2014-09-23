#ifndef COOLFluiD_Numerics_AnalyticalEE_AnalyticEEData_hh
#define COOLFluiD_Numerics_AnalyticalEE_AnalyticEEData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/ErrorEstimatorData.hh"
#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class ConvectiveVarSet; }

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// AnalyticalEECom 's that compose the AnalyticalEE.
/// @author Tiago Quintino
class AnalyticEEData : public Framework::ErrorEstimatorData {

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  AnalyticEEData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~AnalyticEEData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "AnalyticEE";
  }

  Common::SafePtr<Framework::ConvectiveVarSet> getUpdateVarSet() const
  {
    return m_updateVarSet.getPtr();
  }

private: // data

  /// Name of the update variable set
  std::string m_updateVarStr;

  /// Update variable set
  Common::SelfRegistPtr<Framework::ConvectiveVarSet> m_updateVarSet;

}; // end of class AnalyticEEData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for AnalyticalEE
typedef Framework::MethodCommand<AnalyticEEData> AnalyticEECom;

/// Definition of a command provider for AnalyticalEE
typedef Framework::MethodCommand<AnalyticEEData>::PROVIDER AnalyticEEComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_AnalyticEEData_hh
