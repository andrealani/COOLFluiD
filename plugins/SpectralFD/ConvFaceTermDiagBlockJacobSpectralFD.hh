#ifndef COOLFluiD_SpectralFD_ConvFaceTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_ConvFaceTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/ConvFaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the face terms for the
 * spectral finite difference schemes for convection terms, for the diagonal block Jacobian matrices.
 *
 * @author Kris Van Den Abeele
 *
 */
class ConvFaceTermDiagBlockJacobSpectralFD : public ConvFaceTermRHSJacobSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit ConvFaceTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvFaceTermDiagBlockJacobSpectralFD();

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

  /**
 * compute the contribution of the convective face term to both diagonal block Jacobian matrices
   */
//   void computeBothDiagBlockJacobsConvFaceTerm();

  /**
   * compute the contribution of the convective face term to one diagonal block Jacobian matrix
   */
  void computeOneDiagBlockJacobConvFaceTerm(const CFuint side);

  /**
   * add face term contribution to diagonal block Jacobian matrix
   */
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // class ConvFaceTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
