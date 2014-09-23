#ifndef COOLFluiD_SpectralFD_BndFaceTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_BndFaceTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/BndFaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary face terms for the
   * spectral finite difference schemes for both convective and diffusive terms to the diagonal block
   * Jacobian matrix if the full scheme is compact.
   *
   * @author Kris Van Den Abeele
   *
   */
class BndFaceTermDiagBlockJacobSpectralFD : public BndFaceTermRHSJacobSpectralFD {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFDMethodData,BndFaceTermDiagBlockJacobSpectralFD > PROVIDER;

public:

  /**
   * Constructor
   */
  BndFaceTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BndFaceTermDiagBlockJacobSpectralFD();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// Execute on the current TRS
  virtual void executeOnTrs();

  /// compute the contribution of the boundary face term to the Jacobian
  void computeJacobBndFaceTerm();

  /// add boundary face term contribution to diagonal block Jacobian matrix
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // end of class BndFaceTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_BndFaceTermDiagBlockJacobSpectralFD_hh
