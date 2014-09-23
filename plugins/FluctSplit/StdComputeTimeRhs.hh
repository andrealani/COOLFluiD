#ifndef COOLFluiD_Numerics_FluctSplit_StdComputeTimeRhs_hh
#define COOLFluiD_Numerics_FluctSplit_StdComputeTimeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a command that computes the time contribution
/// to the RHS using standard cell center FVM schemes
/// @author Andrea Lani
class FluctSplit_API StdComputeTimeRhs : public FluctuationSplitCom {
public:

  /// Constructor.
  explicit StdComputeTimeRhs(const std::string& name);

  /// Destructor.
  virtual ~StdComputeTimeRhs();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Execute Processing actions
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> _lss;

  /// storage of the States
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// storage of the update coefficients
  Framework::DataSocketSink<
                            CFreal> socket_updateCoeff;

  // the socket to the data handle of arrays of flags specifying if the
  // time jacobian contribution of certain variables in boundary states
  // have to be discarded
  Framework::DataSocketSink<std::vector<bool> > socket_discardTimeJacob;

}; // class StdComputeTimeRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdComputeTimeRhs_hh
