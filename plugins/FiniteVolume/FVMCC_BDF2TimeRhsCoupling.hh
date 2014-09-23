#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_BDF2TimeRhsCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_BDF2TimeRhsCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_PseudoSteadyTimeRhsCoupling.hh"

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
class FVMCC_BDF2TimeRhsCoupling : public FVMCC_PseudoSteadyTimeRhsCoupling {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_BDF2TimeRhsCoupling(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_BDF2TimeRhsCoupling();
  
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
  virtual void computeNumericalTransMatrix(const CFuint iState,
					   const CFuint iLSS);
  
  /**
   * Compute the analytical transformation matrix
   */
  virtual void computeAnalyticalTransMatrix(const CFuint iState,
					    const CFuint iLSS);
  
  
protected: // data

  /// storage of the past time rhs
  Framework::DataSocketSink< CFreal> socket_pastTimeRhs;

}; // class FVMCC_BDF2TimeRhsCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_BDF2TimeRhsCoupling_hh
