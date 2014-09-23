#ifndef COOLFluiD_SpectralFV_DiffFaceTermRHSSpectralFV_hh
#define COOLFluiD_SpectralFV_DiffFaceTermRHSSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFV/BaseFaceTermComputer.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a basic command that computes the face terms for the
 * spectral finite volume schemes for diffusive terms
 *
 * @author Kris Van den Abeele
 *
 */
class DiffFaceTermRHSSpectralFV : public SpectralFVMethodCom {
public:

  /**
   * Constructor.
   */
  explicit DiffFaceTermRHSSpectralFV(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffFaceTermRHSSpectralFV();

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
  void setFaceTermData();

  /// set the gradients of left and right cell
  void setGradients();

  /// add the residual updates to the RHS
  void updateRHS();

  /// add the contributions to the update coefficient
  void addUpdateCoeffContributions();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal > socket_rhs;

  /// socket for updateCoeff
  Framework::DataSocketSink<CFreal > socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// Strategy that computes the face terms for a face
  Common::SafePtr< BaseFaceTermComputer > m_faceTermComputer;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// current face
  Framework::GeometricEntity* m_face;

  /// variable for current neighbouring cells
  std::vector< Framework::GeometricEntity* > m_cells;

  /// CV-CV connectivity through SV face
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_cvCVConnSVFace;

  /// variable for current face orientation
  CFuint m_orient;

  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;

  /// variable for the gradients in the left and right cell
  std::vector< std::vector< std::vector< RealVector >* > > m_grads;

  /// CV face fluxes
  RealVector m_resUpdates;

  /// contributions to the update coefficients
  std::vector< CFreal > m_updateCoefContr;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class DiffFaceTermRHSSpectralFV

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_ComputeConvVolTerms_hh
