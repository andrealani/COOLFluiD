#ifndef COOLFluiD_SpectralFD_DiffVolTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_DiffVolTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/DiffVolTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the volume terms for the
 * spectral finite difference schemes for diffusion terms, for the diagonal block Jacobian matrices.
 *
 * @author Kris Van den Abeele
 *
 */
class DiffVolTermDiagBlockJacobSpectralFD : public DiffVolTermRHSJacobSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit DiffVolTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffVolTermDiagBlockJacobSpectralFD();

  /**
   * Setup private data and data of the aggregated classes
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
   * set the data required to compute the volume term
   */
  void setVolumeTermData();

  /**
   * compute the contribution of the diffusive volume term to the Jacobian
   */
  void computeJacobDiffVolTerm();

  /**
   * add volume term contribution to diagonal block Jacobian matrix
   */
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // class DiffVolTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_DiffVolTermDiagBlockJacobSpectralFD_hh
