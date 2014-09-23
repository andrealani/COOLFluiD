#ifndef COOLFluiD_SpectralFD_VolTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_VolTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/VolTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic command that computes the volume terms for the
 * spectral finite difference schemes for both convective and diffusive terms
 * for the diagonal block Jacobian matrices if the full scheme is compact.
 *
 * @author Kris Van den Abeele
 *
 */
class VolTermDiagBlockJacobSpectralFD : public VolTermRHSJacobSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit VolTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~VolTermDiagBlockJacobSpectralFD();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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

  /// set the required data in the volume and face term computers
  void setVolumeAndFaceTermComputersData();

  /// compute the contribution of the volume term to the Jacobian
  void computeJacobVolTerm();

  /// add volume term contribution to diagonal block Jacobian matrix
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // class VolTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_VolTermDiagBlockJacobSpectralFD_hh
