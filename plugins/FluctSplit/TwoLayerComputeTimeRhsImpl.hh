#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerComputeTimeRhsImpl_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerComputeTimeRhsImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the time contribution
 * to the RHS using standard cell center FVM schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerComputeTimeRhsImpl : public FluctuationSplitCom {
public:

  /**
   * Constructor.
   */
  explicit TwoLayerComputeTimeRhsImpl(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~TwoLayerComputeTimeRhsImpl();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> _lss;

  /// storage of the States
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// storage of the update coefficients (intermediate layer)
  Framework::DataSocketSink<CFreal> socket_interUpdateCoeff;

 }; // class TwoLayerComputeTimeRhsImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerComputeTimeRhsImpl_hh
