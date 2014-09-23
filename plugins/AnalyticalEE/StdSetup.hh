#ifndef COOLFluiD_Numerics_AnalyticalEE_StdSetup_hh
#define COOLFluiD_Numerics_AnalyticalEE_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

#include "AnalyticalEE/AnalyticEEData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to setup the AnalyticalEE
/// @author Tiago Quintino
class StdSetup : public AnalyticEECom
{
public: // methods

  /// Constructor.
  explicit StdSetup(std::string name);

  /// Destructor.
  ~StdSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected: // data

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_StdSetup_hh

