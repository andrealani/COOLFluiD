#ifndef COOLFluiD_FluxReconstructionMethod_BndFaceTermRHSFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_BndFaceTermRHSFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/CompactBndFaceTermComputer.hh"
//#include "FluxReconstructionMethod/BaseBndFaceTermComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary face terms for the
   * spectral finite difference schemes for both convective and diffusive terms to the RHS if the full
   * scheme is compact.
   *
   * @author Kris Van Den Abeele
   *
   */
class BndFaceTermRHSFluxReconstruction : public FluxReconstructionSolverCom {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,BndFaceTermRHSFluxReconstruction > PROVIDER;

public:

  /**
   * Constructor
   */
  BndFaceTermRHSFluxReconstruction(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BndFaceTermRHSFluxReconstruction();

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
  void setFaceTermComputerData();

  /// compute the required data in the face term for the current face
  void computeFaceTermData();

  /// add the residual updates to the RHS
  void addUpdatesToResidual();

  /// add the contributions to the update coefficient
  void addUpdateCoeffContribution();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// Strategy that computes the boundary face terms
  Common::SafePtr< CompactBndFaceTermComputer > m_bndFaceTermComputer; 

  /// the BCStateComputer for this BC
  Common::SafePtr< BCStateComputer > m_bcStateComputer;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current internal cell
  Framework::GeometricEntity* m_intCell;

  /// variable for current face orientation
  CFuint m_orient;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// residual updates
  RealVector m_resUpdates;

  /// diffusive updates to the residuals
  RealVector m_diffResUpdates;

  /// contribution to the update coefficient in the neighbouring cell
  CFreal m_updateCoeffContr;

  /// diffusive contribution to the update coefficient in the neighbouring cell
  CFreal m_diffUpdateCoeffContr;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // end of class BndFaceTermRHSFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BndFaceTermRHSFluxReconstruction_hh
