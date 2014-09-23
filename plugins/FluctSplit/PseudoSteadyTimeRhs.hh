#ifndef COOLFluiD_Numerics_FluctSplit_PseudoSteadyTimeRhs_hh
#define COOLFluiD_Numerics_FluctSplit_PseudoSteadyTimeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit//StdComputeTimeRhs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a command that computes the pseudo steady RHS
/// using RD schemes
/// @author Andrea Lani
class FluctSplit_API PseudoSteadyTimeRhs : public StdComputeTimeRhs {
public:

  /// Constructor.
  explicit PseudoSteadyTimeRhs(const std::string& name);

  /// Destructor.
  ~PseudoSteadyTimeRhs();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Execute Processing actions
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;

  /// storage of the past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;

  /// storage of the rhs
  Framework::DataSocketSink<
                            CFreal> socket_rhs;

  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;

  /// temporary state
  RealVector _tempState;

  /// temporary perturbed state
  RealVector _tempPertState;

  // accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> _acc;

}; // class PseudoSteadyTimeRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PseudoSteadyTimeRhs_hh
