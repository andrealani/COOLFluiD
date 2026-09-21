// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the diffusive part of the system using a FluxReconstruction solver for an implicit scheme
/// @author Alexander Papen
/// @author Ray Vandenhoeck
/// @author Rayan Dhib
class DiffRHSJacobFluxReconstruction : public DiffRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit DiffRHSJacobFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~DiffRHSJacobFluxReconstruction() {}

  /// Execute processing actions
  virtual void execute();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();
    
protected: //functions

  /**
   * set the local indexes of the other faces (not the current boundary face)
   * @pre m_faces is set
   */
  void setOtherFacesLocalIdxs();
  
  /**
   * compute the unperturbed cell diffusive residuals
   * @pre m_faceTermComputers->computeDiffFaceTermAndUpdateCoefContributions
   * @pre setCellsData()
   */
  virtual void computeUnpertCellDiffResiduals(const CFuint side);
  
  /**
   * Add to m_cellGrads the change of the all-face gradients of both cells of
   * the current face when the state m_pertVar of solution point m_pertSol on
   * side is perturbed: the volume term and the corrections of every face of
   * the perturbed cell, the correction of the current face also in the other
   * cell.
   */
  virtual void computePerturbedGradientsAnalytical(const CFuint side);

  /**
   * Add to m_cellGrads the change of the correction of the current face in both
   * cells: the jump changes by -0.5*dg^D_f in the perturbed cell and by
   * +0.5*dg^D_f in the other cell, with dg^D_f the change of the gradient
   * variables extrapolated to the flux points of the face.
   * @pre addPerturbedVolumeGradient(side)
   */
  void addPerturbedCurrentFaceGradient(const CFuint side);

  /**
   * Compute the boundary value of the gradient variables at the flux points of
   * the current boundary face: the rule of the boundary condition for the
   * physical gradient variables, 0.5*(g^D_f + g_AV(U_ghost)) for the artificial
   * viscosity variables.
   * @pre the ghost states are in m_flxPntGhostSol
   * @param gradVarsFlxPnt gradient variables extrapolated to the flux points
   * @param bndGradVars output boundary value of the gradient variables
   * @param artificialViscosity true for the artificial viscosity variables
   */
  void computePertBndGradVars(Common::SafePtr< BCStateComputer > bc,
                              const std::vector< RealVector* >& gradVarsFlxPnt,
                              std::vector< RealVector* >& bndGradVars,
                              const bool artificialViscosity);

  /**
   * Add to m_cellGrads[side] the change of the volume term of the gradient of
   * the perturbed cell, and set m_pertGradVarsChange to the change of its
   * gradient variables at the perturbed solution point.
   * @param artificialViscosity true for the artificial viscosity variables
   */
  void addPerturbedVolumeGradient(const CFuint side, const bool artificialViscosity = false);

  /**
   * Add to m_cellGrads[side] the change of the corrections of the given faces
   * of the perturbed cell: the jump changes by -0.5*dg^D_f on an interior face
   * and by dg_b - dg^D_f on a boundary face, with dg^D_f the change of the
   * gradient variables extrapolated to the flux points and dg_b the change of
   * their boundary value.
   * @pre addPerturbedVolumeGradient(side)
   * @param artificialViscosity true for the artificial viscosity variables
   */
  void addPerturbedFaceLiftings(const CFuint side, const std::vector< CFuint >& faceLocalIdxs, const bool artificialViscosity = false);

  /**
   * Volume residual and its Jacobian block of a cell whose faces are all
   * boundary faces, which the face loop never visits.
   * @pre m_cells[LEFT] and m_states[LEFT] hold the cell, built with its faces
   */
  void computeCellWithoutInnerFace(const CFuint cellID, const bool artificialViscosity = false);

  /**
   * Prepare the command for the perturbations of a cell whose faces are all
   * boundary faces, after its residual is computed.
   */
  virtual void prepareIsolatedCellJacobian()
  {
  }

  /**
   * compute the perturbed cell diffusive residuals for one cell
   * @pre m_faceTermComputer->computeDiffFaceTerm
   * @pre backupAndReconstructOtherFacesAndCellPhysVars(
   * @pre reconstructOtherFacesAndCellGradients()
   */
  virtual void computePertCellDiffResiduals(const CFuint side);

  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  void computeOneJacobDiffFaceTerm(const CFuint side);
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);

  /**
   * Extrapolate the gradient variables of the perturbed cell to the flux points
   * of one face, with the perturbation (m_pertGradVarsFlxPnt) and without it
   * (m_gradVarsFlxPntBefore).
   * @param flxPntConn cell flux point index of every flux point of the face
   */
  void extrapolateGradVarsToFaceFlxPnts(const std::vector< CFuint >& flxPntConn);

  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);

  /**
   * compute the terms for the gradient computation for one side of a face,
   * from the states currently held in m_cellStatesFlxPnt[side]. The face jump
   * terms take the gradient variables of the extrapolated face state, so their
   * linearisation uses the change of the gradient variables at the face, not
   * the change at the solution point times the basis value. The two differ
   * when the gradient variables are not affine in the state (mass fractions,
   * p/rho).
   */
  virtual void computeFlxPntGradTerm(const CFuint side, RealMatrix& gradTerm);

  /**
   * compute the data needed for the computation of the perturbed gradients
   */
  void computePertGradData(const CFuint side);
  
  /**
   * store values that will be overwritten
   */
  void storeBackups();
  
  /**
   * restore values that were overwritten
   */
  void restoreFromBackups();
  
  /**
   *  add the residual updates to the RHS
   */
  void updateRHSUnpertCell(const CFuint side);
  
  /// compute the divergence of the discontinuous flux (-divFD+divhFD) of a neighbor cell
  virtual void computeDivDiscontFlxNeighb(RealVector& residuals, const CFuint side);
  
  
protected: //data
  
  /// builder of cells
  std::vector< Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > > m_cellBuilders;
  
  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// single cell accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_accSC;

  /// variable for faces
  std::vector< const std::vector< Framework::GeometricEntity* >* > m_faces;

  /// perturbed updates to the residuals
  std::vector< RealVector > m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// unperturbed diffusive residuals
  std::vector< RealVector > m_unpertCellDiffRes;

  /// perturbed diffusive residuals
  RealVector m_pertCellDiffRes;

  /// derivative diffusive residuals
  RealVector m_derivCellDiffRes;

  /// cell local indexes of the other faces (not the face itself)
  std::vector< std::vector< CFuint > > m_otherFaceLocalIdxs;

  /// pointer to booleans telling whether a face is on the boundary
  std::vector< Common::SafePtr< std::vector< bool > > > m_isFaceOnBoundary;

  /// pointer to neighbouring cell side vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_nghbrCellSide;

  /// pointer to current cell side vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_currCellSide;

  /// pointer to orientation vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_faceOrients;

  /// pointer to BC index vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_faceBCIdx;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;
  
  /// ghost flux point solutions
  std::vector< Framework::State* > m_flxPntGhostSol;
  
  /// Divergence of the continuous flux at the solution points of the left neighbour
  std::vector< RealVector> m_divContFlxL;
  
  /// Divergence of the continuous flux at the solution points of the right neighbour
  std::vector< RealVector> m_divContFlxR;
  
  /// unperturbed updates to the residuals
  std::vector< RealVector > m_resUpdates;
  
  /// backup of the gradients in the neighbouring cell
  std::vector< std::vector< std::vector< RealVector > > > m_cellGradsBackUp;
  
  /// Perturbed divergence of the continuous flux at the solution points of the neighbours
  std::vector< std::vector< RealVector> > m_pertDivContFlx;
  
  /// perturbed corrections
  std::vector< RealVector> m_pertCorrections;
  
  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell;

  /// pointer to neighbouring cell side vector
  Common::SafePtr< std::vector< CFuint > > m_nghbrCellSideCell;

  /// pointer to current cell side vector
  Common::SafePtr< std::vector< CFuint > > m_currCellSideCell;

  /// pointer to orientation vector
  Common::SafePtr< std::vector< CFuint > > m_faceOrientsCell;

  /// pointer to BC index vector
  Common::SafePtr< std::vector< CFuint > > m_faceBCIdxCell;
  
  /// the ghost gradients in the flux points
  std::vector< std::vector< RealVector* > > m_flxPntGhostGrads;
  
  /// the current flux pnt number
  CFuint m_currFlx;
  
  /// list of the vectors to which to calculate the derivative
  std::vector< std::vector< CFuint > > m_dimList;
  
  /// gradient variables of the perturbed cell, with the perturbation
  RealMatrix m_pertGradVarsSolPnts;
  
  /// gradient variables of the perturbed cell, without the perturbation
  RealMatrix m_gradVarsSolPntsBefore;

  /// storage of the gradient variables extrapolated to the flux points of a face, perturbed
  std::vector< RealVector > m_pertGradVarsFlxPntStore;
  
  /// storage of the gradient variables extrapolated to the flux points of a face, unperturbed
  std::vector< RealVector > m_gradVarsFlxPntBeforeStore;
  
  /// storage of the boundary value of the gradient variables, perturbed
  std::vector< RealVector > m_pertBndGradVarsStore;
  
  /// storage of the boundary value of the gradient variables, unperturbed
  std::vector< RealVector > m_bndGradVarsBeforeStore;
  
  /// gradient variables extrapolated to the flux points of a face, perturbed
  std::vector< RealVector* > m_pertGradVarsFlxPnt;
  
  /// gradient variables extrapolated to the flux points of a face, unperturbed
  std::vector< RealVector* > m_gradVarsFlxPntBefore;
  
  /// boundary value of the gradient variables at the flux points, perturbed
  std::vector< RealVector* > m_pertBndGradVars;
  
  /// boundary value of the gradient variables at the flux points, unperturbed
  std::vector< RealVector* > m_bndGradVarsBefore;
  
  /// vector to temporarily store a correction projected on a normal
  RealVector m_projectedCorrL;
  
  /// vector to temporarily store a correction projected on a normal
  RealVector m_projectedCorrR;
  
  /// perturbed side
  CFuint m_pertSide;
  
  /// perturbed sol pnt
  CFuint m_pertSol;
  
  /// perturbed variable
  CFuint m_pertVar;
  
  /// perturbation value
  RealVector m_pertGradVarsChange;

  /// flags for each cell to tell whether its inner -divFD has been computed 
  std::vector< bool > m_cellFlags;
  
  /// unperturbed diffusive residuals for all cells
  std::vector< RealVector > m_unpertAllCellDiffRes;
  
  /// bools telling whether solution points are affected by perturbation (for both neighbor cells)
  std::vector< std::vector < bool > > m_affectedSolPnts;
  
  /// Continuous flux at the solution points backup for both neighbor cells
  std::vector< std::vector< std::vector< RealVector > > > m_contFlxBackup;
  
  /// Continuous flux at the solution points for both neighbor cells
  std::vector< std::vector< std::vector< RealVector > > > m_contFlxNeighb;


  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstruction_hh

