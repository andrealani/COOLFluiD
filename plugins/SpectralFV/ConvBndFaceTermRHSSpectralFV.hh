#ifndef COOLFluiD_SpectralFV_ConvBndFaceTermRHSSpectralFV_hh
#define COOLFluiD_SpectralFV_ConvBndFaceTermRHSSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/SpectralFVElementData.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary face terms for the
   * spectral finite volume schemes for convective terms to the RHS
   *
   * @author Kris Van Den Abeele
   *
   */
class ConvBndFaceTermRHSSpectralFV : public SpectralFVMethodCom {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFVMethodData,ConvBndFaceTermRHSSpectralFV > PROVIDER;

public:

  /**
   * Constructor
   */
  ConvBndFaceTermRHSSpectralFV(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvBndFaceTermRHSSpectralFV();

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

  /**
   * set the BC state computer
   */
  void setBcStateComputer(Common::SafePtr< BCStateComputer > bcStateComputer)
  {
    m_bcStateComputer = bcStateComputer;
  }

protected: // functions

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();

  /// set face term data for current element type
  void setFaceTermData();

  /// add the residual updates to the RHS
  void updateRHS();

  /// add the updates to the wave speed
  void updateWaveSpeed();

  /// compute the face term contribution to the gradients
  void computeGradientFaceTerm();

  /// add updates to gradients
  void addGradBCTerms();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// Strategy that computes the boundary face terms
  Common::SafePtr< BaseBndFaceTermComputer > m_bndFaceTermComputer;

  /// the BCStateComputer for this BC
  Common::SafePtr< BCStateComputer > m_bcStateComputer;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current cell
  Framework::GeometricEntity* m_intCell;

  /// boundary face CV connectivity on SV face
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_svFaceCVConn;

  /// variable for current face orientation
  CFuint m_orient;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// CV face fluxes
  RealVector m_resUpdates;

  /// update for the wave speed in the neighbouring cell
  CFreal m_waveSpeedUpd;

  /// CV face gradient updates
  std::vector< std::vector< RealVector > > m_gradUpdates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // end of class ConvBndFaceTermRHSSpectralFV

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_ConvBndFaceTermRHSSpectralFV_hh
