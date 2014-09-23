#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsTridiag_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsTridiag_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_PseudoSteadyTimeRhs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the pseudo steady RHS
 * using standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_PseudoSteadyTimeRhsTridiag : public FVMCC_PseudoSteadyTimeRhs {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_PseudoSteadyTimeRhsTridiag(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_PseudoSteadyTimeRhsTridiag();
 
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /**
   * Compute the transformation matrix numerically
   */
  virtual void computeNumericalTransMatrix(const CFuint iState);

  /**
   * Compute the analytical transformation matrix
   */
  virtual void computeAnalyticalTransMatrix(const CFuint iState);

private:

  /// storage of diagonal block matrices  
  Framework::DataSocketSink<CFreal> socket_diagMatrices; 

  /// storage of the local updatable IDs or -1 (ghost) for all local states
  Framework::DataSocketSink<CFint> socket_upLocalIDsAll;

  /// matrix iterator
  RealMatrix _matIter;

}; // class FVMCC_PseudoSteadyTimeRhsTridiag

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsTridiag_hh
