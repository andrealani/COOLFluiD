#ifndef COOLFluiD_SpectralFD_FaceTermRHSSpectralFD_hh
#define COOLFluiD_SpectralFD_FaceTermRHSSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/CompactFaceTermComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a basic command that computes the contribution of face terms for the
 * spectral finite difference schemes for both convective and diffusive terms to the RHS if the full
 * scheme is compact.
 *
 * @author Kris Van den Abeele
 *
 */
class FaceTermRHSSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit FaceTermRHSSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FaceTermRHSSpectralFD();

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

  /// set face term data
  void setFaceTermComputerData();

  /// compute the required data in the face term for the current face
  void computeFaceTermData();

  /// add the residual updates to the RHS
  void addUpdatesToResidual();

  /// add the contributions to the update coefficient
  void addUpdateCoeffContributions();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal > socket_rhs;

  /// socket for updateCoeff
  Framework::DataSocketSink<CFreal > socket_updateCoeff;

  /// Strategy that computes the face terms for a face
  Common::SafePtr< CompactFaceTermComputer > m_faceTermComputer;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// current face
  Framework::GeometricEntity* m_face;

  /// variable for current neighbouring cells
  std::vector< Framework::GeometricEntity* > m_cells;

  /// variable for current face orientation
  CFuint m_orient;

  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;

  /// residual updates
  std::vector< RealVector > m_resUpdates;

  /// diffusive residual updates
  std::vector< RealVector > m_diffResUpdates;

  /// contributions to the update coefficients
  std::vector< CFreal > m_updateCoefContr;

  /// diffusive contributions to the update coefficients
  std::vector< CFreal > m_diffUpdateCoefContr;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class FaceTermRHSSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
