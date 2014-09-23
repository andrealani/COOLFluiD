#ifndef COOLFluiD_SpectralFD_DiffBndFaceTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_DiffBndFaceTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/DiffBndFaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////


/**
 * This class represents a command that computes contribution of the boundary face terms for the
 * spectral finite volume schemes for diffusion terms, for the diagonal block Jacobian matrices.
 *
 * @author Kris Van Den Abeele
 *
 */
class DiffBndFaceTermDiagBlockJacobSpectralFD : public DiffBndFaceTermRHSJacobSpectralFD {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFDMethodData,DiffBndFaceTermDiagBlockJacobSpectralFD > PROVIDER;

public:

  /**
   * Constructor
   */
  DiffBndFaceTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndFaceTermDiagBlockJacobSpectralFD();

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
   * compute the contribution of the diffusive boundary face term to the Jacobian
   */
  void computeJacobDiffBndFaceTerm();

  /**
   * add boundary face term contribution to diagonal block Jacobian matrix
   */
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // end of class DiffBndFaceTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_DiffBndFaceTermDiagBlockJacobSpectralFD_hh
