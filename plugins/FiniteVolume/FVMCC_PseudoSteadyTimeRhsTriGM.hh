#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsTriGM_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsTriGM_hh

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
 * @author Jiri Simonek
 *
 */
class FVMCC_PseudoSteadyTimeRhsTriGM : public FVMCC_PseudoSteadyTimeRhs {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_PseudoSteadyTimeRhsTriGM(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_PseudoSteadyTimeRhsTriGM();

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

}; // class FVMCC_PseudoSteadyTimeRhsTriGM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsTriGM_hh
