#ifndef COOLFluiD_SpectralFD_LUSGSBDF3Prepare_hh
#define COOLFluiD_SpectralFD_LUSGSBDF3Prepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/LUSGSPrepare.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to prepare
 * the computation for the LU-SGS algorithm with BDF3 scheme
 * and variable time step.
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class LUSGSBDF3Prepare : public LUSGSPrepare {

public: // functions

  /**
   * Constructor.
   */
  explicit LUSGSBDF3Prepare(const std::string& name);

  /**
   * Destructor.
   */
  ~LUSGSBDF3Prepare();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:

  /// pointer to the 3 steps time marching scheme parameters (BDF3 with variable time step)
  Common::SafePtr< RealVector > m_3StepsTMSparams;

}; // class Prepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_LUSGSBDF3Prepare_hh
