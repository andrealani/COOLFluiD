#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ALEBDF2TimeRhsCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ALEBDF2TimeRhsCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ALETimeRhsCoupling.hh"

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
 * @author Thomas Wuilbaut
 *
 */
class FVMCC_ALEBDF2TimeRhsCoupling : public FVMCC_ALETimeRhsCoupling {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ALEBDF2TimeRhsCoupling(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_ALEBDF2TimeRhsCoupling();

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
  void computeNumericalTransMatrix(const CFuint iState,
			           const CFuint iLSS);

  /**
   * Compute the analytical transformation matrix
   */
  void computeAnalyticalTransMatrix(const CFuint iState,
			            const CFuint iLSS);


protected: // data

  /// storage of the past time rhs
  Framework::DataSocketSink< CFreal> socket_pastTimeRhs;

}; // class FVMCC_ALEBDF2TimeRhsCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ALEBDF2TimeRhsCoupling_hh
