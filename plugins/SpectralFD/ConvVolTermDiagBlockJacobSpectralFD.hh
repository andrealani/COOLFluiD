#ifndef COOLFluiD_SpectralFD_ConvVolTermDiagBlockJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_ConvVolTermDiagBlockJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/ConvVolTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the volume terms for the
 * spectral finite difference schemes for convection terms, for the diagonal block Jacobian matrices.
 *
 * @author Kris Van den Abeele
 *
 */
class ConvVolTermDiagBlockJacobSpectralFD : public ConvVolTermRHSJacobSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit ConvVolTermDiagBlockJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvVolTermDiagBlockJacobSpectralFD();

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

  /**
   * set the data required to compute the volume term
   */

  void setVolumeTermData();

  /**
   * compute the contribution of the convective volume term to the diagonal block Jacobian matrixx
   */
  void computeDiagBlockJacobConvVolTerm();

  /**
   * add volume term contribution to diagonal block Jacobian matrix
   */
  void addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx);

protected: // data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// current diagonal block Jacobian matrix
RealMatrix* m_currDiagMatrix;

}; // class ConvVolTermDiagBlockJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ConvVolTermDiagBlockJacobSpectralFD_hh
