#ifndef COOLFluiD_SpectralFV_DiffBndFaceTermRHSSpectralFV_hh
#define COOLFluiD_SpectralFV_DiffBndFaceTermRHSSpectralFV_hh

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
   * spectral finite volume schemes for diffusive terms to the RHS
   *
   * @author Kris Van Den Abeele
   *
   */
class DiffBndFaceTermRHSSpectralFV : public SpectralFVMethodCom {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFVMethodData,DiffBndFaceTermRHSSpectralFV > PROVIDER;

public:

  /**
   * Constructor
   */
  DiffBndFaceTermRHSSpectralFV(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndFaceTermRHSSpectralFV();

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

  /// set the gradients of left and right cell
  void setGradients();

  /// add the residual updates to the RHS
  void updateRHS();

  /// add the contributions to the update coefficient
  void addUpdateCoeffContribution();

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

  /// variable for current internal cell
  Framework::GeometricEntity* m_intCell;

  /// boundary face CV connectivity on SV face
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_svFaceCVConn;

  /// variable for current face orientation
  CFuint m_orient;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// the gradients in the neighbouring cell
  std::vector< std::vector< RealVector >* > m_cellGrads;

  /// CV face fluxes
  RealVector m_resUpdates;

  /// contribution to the update coefficient in the neighbouring cell
  CFreal m_updateCoeffContr;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // end of class DiffBndFaceTermRHSSpectralFV

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_DiffBndFaceTermRHSSpectralFV_hh
