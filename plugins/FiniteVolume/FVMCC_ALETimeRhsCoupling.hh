#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ALETimeRhsCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ALETimeRhsCoupling_hh

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
class FVMCC_ALETimeRhsCoupling : public FVMCC_PseudoSteadyTimeRhsCoupling {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ALETimeRhsCoupling(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ALETimeRhsCoupling();

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

  /// storage of past volumes
  Framework::DataSocketSink< CFreal> socket_pastVolumes;

  /// past part of the diagonal value in the block to be inserted in the LSS matrix
  CFreal _pastDiagValue;

  CFreal _volume;
  CFreal _dt;

}; // class FVMCC_ALETimeRhsCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ALETimeRhsCoupling_hh
