#ifndef COOLFluiD_SpectralFD_DiffFaceTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_DiffFaceTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/DiffFaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the face terms for the
 * spectral finite volume schemes for diffusion terms, for the diagonal block Jacobian matrices.
 *
 * @author Kris Van den Abeele
 *
 */
class DiffFaceTermDiagBlockJacobSpectralFD : public DiffFaceTermRHSJacobSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit DiffFaceTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffFaceTermDiagBlockJacobSpectralFD();

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
   * compute the contribution of the diffusive face term to both Jacobians
   */
//   void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  void computeOneJacobDiffFaceTerm(const CFuint side);

  /**
   * add face term contribution to diagonal block Jacobian matrix
   */
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrices
RealMatrix* m_currDiagMatrix;

}; // class DiffFaceTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_DiffFaceTermDiagBlockJacobSpectralFD_hh
