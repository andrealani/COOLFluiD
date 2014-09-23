#ifndef COOLFluiD_SpectralFD_ConvBndFaceTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_ConvBndFaceTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/ConvBndFaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////


/**
  * This class represents a command that computes contribution of the boundary face terms for the
  * spectral finite volume schemes for convection terms, for the diagonal block Jacobian matrices.
  *
  * @author Kris Van Den Abeele
  *
  */
class ConvBndFaceTermDiagBlockJacobSpectralFD : public ConvBndFaceTermRHSJacobSpectralFD {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFDMethodData,ConvBndFaceTermDiagBlockJacobSpectralFD > PROVIDER;

public:

  /**
   * Constructor
   */
  ConvBndFaceTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvBndFaceTermDiagBlockJacobSpectralFD();

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

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();

  /**
   * compute the contribution of the convective boundary face term to the diagonal block Jacobian matrix
   */
  void computeDiagBlockJacobConvBndFaceTerm();

  /**
   * add boundary face term contribution to diagonal block Jacobian matrix
   */
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // end of class ConvBndFaceTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ConvBndFaceTermDiagBlockJacobSpectralFD_hh
