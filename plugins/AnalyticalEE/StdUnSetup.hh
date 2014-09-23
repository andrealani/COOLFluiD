#ifndef COOLFluiD_Numerics_AnalyticalEE_StdUnSetup_hh
#define COOLFluiD_Numerics_AnalyticalEE_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "AnalyticEEData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Tiago Quintino
class StdUnSetup : public AnalyticEECom {
public:

  /// Constructor.
  explicit StdUnSetup(std::string name);

  /// Destructor.
  ~StdUnSetup();

  /// Execute Processing actions
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_StdUnSetup_hh

