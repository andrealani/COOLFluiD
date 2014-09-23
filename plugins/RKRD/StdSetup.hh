#ifndef COOLFluiD_Numerics_RKRD_StdSetup_hh
#define COOLFluiD_Numerics_RKRD_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKRDData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "MathTools/RealVector.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to setup the RKRD Method
class StdSetup : public RKRDCom {
public:

  /// Constructor.
  explicit StdSetup(const std::string& name);

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  virtual void execute();

protected:

  /// socket for Rhs
  Framework::DataSocketSource<CFreal> socket_rhs;
  /// socket for updateCoeff
  Framework::DataSocketSource<CFreal> socket_updateCoeff;
  /// socket for states at different k-steps
  Framework::DataSocketSource<RealMatrix> socket_kstates;
  /// socket with median cell areas
  Framework::DataSocketSource<CFreal> socket_median_areas;
  /// socket to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RKRD_StdSetup_hh

