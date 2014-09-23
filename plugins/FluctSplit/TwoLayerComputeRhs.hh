#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerComputeRhs_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerComputeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeRHS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class State; }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a command that computes the RhsAndJacob using
/// standard RD schemes and Petsc library
/// @author Thomas Wuilbaut
class TwoLayerComputeRhs : public ComputeRHS {
public:

  /// Constructor.
  explicit TwoLayerComputeRhs(const std::string& name);

  /// Destructor.
  virtual ~TwoLayerComputeRhs();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Cleans the rhs setting it to zero
  void cleanRHS()
  {
    Framework::DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
    Framework::DataHandle< CFreal> interRhs = socket_interRhs.getDataHandle();
    rhs = 0.0;
    interRhs = 0.0;
  }

protected: // functions

  /// Execute the command on the current TRS
  virtual void executeOnTrs();

protected: // variables

  /// temporary storage
  std::vector<RealVector> _temp;

  /// handle for the interRhs
  Framework::DataSocketSink<CFreal>   socket_interRhs;

  /// handle for the interRhs
  Framework::DataSocketSink<CFreal>   socket_interUpdateCoeff;

}; // class TwoLayerComputeRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerComputeRhs_hh
