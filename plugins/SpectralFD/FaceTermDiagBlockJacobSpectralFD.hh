#ifndef COOLFluiD_SpectralFD_FaceTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_FaceTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/FaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a basic command that computes the contribution of face terms for the
 * spectral finite difference schemes for both convective and diffusive terms to the diagonal block Jacobian
 * if the full scheme is compact.
 *
 * @author Kris Van den Abeele
 *
 */
class FaceTermDiagBlockJacobSpectralFD : public FaceTermRHSJacobSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit FaceTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FaceTermDiagBlockJacobSpectralFD();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// compute the contribution of the face term to one Jacobians
  void computeOneJacobFaceTerm(const CFuint side);

  /// add face term contribution to diagonal block Jacobian matrix
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // class FaceTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
