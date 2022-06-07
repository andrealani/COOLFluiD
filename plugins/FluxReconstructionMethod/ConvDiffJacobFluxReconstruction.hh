// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the convective, diffusive 
/// and artificial viscosity part of the system using a FluxReconstruction solver for an implicit scheme
/// @author Ray Vandenhoeck
class ConvDiffJacobFluxReconstruction : public DiffRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvDiffJacobFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~ConvDiffJacobFluxReconstruction() {}

  /// Execute processing actions
  virtual void execute();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
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
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
    
protected: //functions
    
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * Set the data for the current face necessary to calculate FI
   */
  virtual void setFaceData(CFuint faceID);
  
  /**
   * Compute the left and right states and gradients in the flx pnts
   */
  virtual void computeFlxPntStatesAndGrads();
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * compute the unperturbed cell diffusive residuals
   * @pre m_faceTermComputers->computeDiffFaceTermAndUpdateCoefContributions
   * @pre setCellsData()
   */
  virtual void computeUnpertCellDiffResiduals(const CFuint side);

  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  void computeOneJacobDiffFaceTerm(const CFuint side);
  
  /**
   * compute the artificial diffusive flux
   */
  virtual void computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux);
  
  /// compute the volume term contribution to the gradients
  virtual void computeGradients();
  
  /// compute the face correction to the corrected gradients
  virtual void computeGradientFaceCorrections();
  
  /// set face data for gradient computation
  void setFaceDataForGradients(const CFuint faceID);
  
  /// initialize the data needed for the jacobian
  virtual void initJacobianComputation();
  
  /// compute the cell flux jacobian numerically
  virtual void computeCellFluxJacobianNum(const CFreal resFactor);
  
  /// compute the Riemann flux jacobian numerically
  virtual void computeRiemannFluxJacobianNum(const CFreal resFactor);
  
  /// compute the flux to gradient jacobian numerically
  virtual void computeFluxToGradJacobianNum(const CFreal resFactor);
  
  /// compute the gradient to state jacobian analytically
  virtual void computeGradToStateJacobianAna();
  
  /// compute the Riemann flux to gradient jacobian numerically
  virtual void computeRiemannFluxToGradJacobianNum(const CFreal resFactor);
  
  /// compute the gradient variables to state jacobians numerically
  virtual void computeGradVarsToStateJacobianNum();
  
protected: //data
    
  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;
  
  /// current element index
  CFuint m_elemIdx;
  
  /// variable for faces
  const std::vector< Framework::GeometricEntity* >* m_facesCell;
  
  /// bool telling whether the jacobian is being computed
  bool m_jacob;
  
  /// storage for the normals in the solution points
  Framework::DataSocketSink< CFreal > socket_solPntNormals;
  
  /// storage for the normals in the solution points
  Framework::DataSocketSink< CFreal > socket_flxPntNormals;
  
  /// storage for the cell volumes
  Framework::DataSocketSink< CFreal > socket_cellVolumes;
  
  /// storage of the geometrical jacobians in the states
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// vector to store sol pnt values temporarily
  RealVector m_tempSolPntVec;
  
  /// vector to store sol pnt values temporarily
  RealVector m_tempSolPntVec2;
  
  /// stores the flux jacobian for each side, in each sol pnt, for each variable, for each direction
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_fluxJacobian;
  
  /// stores the Riemann flux jacobian for each side, in each face flx pnt, for each variable
  std::vector< std::vector< std::vector< RealVector > > > m_riemannFluxJacobian;
  
  /// convective Riemann flux
  std::vector < RealVector > m_flxPntRiemannFluxDiff;
  
  /// perturbed convective Riemann Flux
  std::vector < RealVector > m_flxPntRiemannFluxPert;
  
  /// temporary storage for a flux
  RealVector m_tempFlux;
  
  /// stores the flux jacobian to the gradients for each side, in each sol pnt, for each variable, for each gradient direction for each flux direction
  std::vector< std::vector< std::vector< std::vector< std::vector< RealVector > > > > > m_gradientFluxJacobian;
  
  /// stores the gradient variables jacobian to the states for each side, in each sol pnt, for each depending state variable, for each grad vars variable
  std::vector< std::vector< std::vector< RealVector > > > m_gradVarsToStateJacobian;
  
  /// stores the gradient jacobian to the states for each side, in each sol pnt, for each depending side for each depending sol pnt, for each gradient direction
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_gradientStateJacobian;
  
  /// stores the Riemann flux jacobian to the gradients for each face flx pnt, for each variable, for each gradient direction
  std::vector< std::vector< std::vector< RealVector > > > m_riemannFluxGradJacobian;
  
  /// Continuous diffusive flux at the solution points backup for both neighbor cells
  std::vector< std::vector< std::vector< RealVector > > > m_contFlxBackupDiff;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // class ConvDiffJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstruction_hh
