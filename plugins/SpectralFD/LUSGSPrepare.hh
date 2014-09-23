#ifndef COOLFluiD_SpectralFD_LUSGSPrepare_hh
#define COOLFluiD_SpectralFD_LUSGSPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to prepare
 * the computation for the LU-SGS algorithm by clearing the diagonal block Jacobian matrices.
 *
 * @author Kris Van den Abeele
 */
class LUSGSPrepare : public SpectralFDMethodCom {

public: // functions

  /**
   * Constructor.
   */
  explicit LUSGSPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~LUSGSPrepare();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data

  /// socket for diagonal block Jacobian matrices
  Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink< CFreal > socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

}; // class Prepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_LUSGSPrepare_hh
