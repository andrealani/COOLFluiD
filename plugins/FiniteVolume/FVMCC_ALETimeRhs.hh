#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ALETimeRhs_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ALETimeRhs_hh

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
class FVMCC_ALETimeRhs : public FVMCC_PseudoSteadyTimeRhs {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ALETimeRhs(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ALETimeRhs();

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
  virtual void computeNumericalTransMatrix(const CFuint iState);

  /**
   * Compute the analytical transformation matrix
   */
  virtual void computeAnalyticalTransMatrix(const CFuint iState);

protected: // data

  /// storage of past volumes
  Framework::DataSocketSink< CFreal> socket_pastVolumes;

  /// past part of the diagonal value in the block to be inserted in the LSS matrix
  CFreal _pastDiagValue;

  CFreal _volume;
  CFreal _dt;

}; // class FVMCC_ALETimeRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ALETimeRhs_hh
