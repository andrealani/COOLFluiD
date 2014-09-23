#ifndef COOLFluiD_Numerics_FluctSplit_PseudoSteadyTimeAnalytMatRhs_hh
#define COOLFluiD_Numerics_FluctSplit_PseudoSteadyTimeAnalytMatRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/StdComputeTimeRhs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace MathTools {
    class MatrixInverter;
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a command that computes the pseudo steady RHS
/// using RD schemes
/// @author Andrea Lani
class FluctSplit_API PseudoSteadyTimeAnalytMatRhs : public StdComputeTimeRhs {
public:

  /// Constructor.
  explicit PseudoSteadyTimeAnalytMatRhs(const std::string& name);

  /// Destructor.
  ~PseudoSteadyTimeAnalytMatRhs();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Execute Processing actions
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       _inverter;

  /// socket to the rhs
  Framework::DataSocketSink<
                            CFreal> socket_rhs;

  /// socket to the pastStates
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;

  // accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> _acc;

}; // class PseudoSteadyTimeAnalytMatRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PseudoSteadyTimeAnalytMatRhs_hh
