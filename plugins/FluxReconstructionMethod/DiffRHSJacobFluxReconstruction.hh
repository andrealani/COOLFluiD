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
   * recompute the cell gradients from the current cell solutions,
   * after perturbation
   */
  virtual void computePerturbedGradientsAnalytical(const CFuint side);

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
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);
  
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

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< std::vector< Framework::State* >* > > > m_faceNghbrStates;

  /// vector containing pointers to the left and right gradients with respect to a face
  std::vector< std::vector< std::vector< std::vector< std::vector< RealVector >* > > > > m_faceNghbrGrads;

  /// perturbed updates to the residuals
  std::vector< RealVector > m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// perturbed left and right cell gradients
  std::vector< std::vector< std::vector< RealVector >* > > m_pertGrads;

  /// left and right cell gradients minus current face term
  std::vector< std::vector< std::vector< RealVector > > > m_cellGradsMinusFaceTerm;

  /// left and right cell gradients minus other face terms
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_cellGradsMinusOtherFaceTerm;

  /// unperturbed diffusive residuals
  std::vector< RealVector > m_unpertCellDiffRes;

  /// perturbed diffusive residuals
  RealVector m_pertCellDiffRes;

  /// derivative diffusive residuals
  RealVector m_derivCellDiffRes;

  /// Jacobian determinants
  std::vector< std::valarray<CFreal> > m_solJacobDet;

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
  
  /// transformed states in a left cell for gradient computation
  RealMatrix m_gradTermL;
  
  /// transformed states in a right cell for gradient computation
  RealMatrix m_gradTermR;
  
  /// temp transformed states in a right cell for gradient computation
  RealMatrix m_gradTermTemp;
  
  /// transformed states in a cell for gradient computation
  RealMatrix m_gradTerm;
  
  /// transformed states in a cell for gradient computation before perturbation
  RealMatrix m_gradTermBefore;
  
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
  
  /// the corrected gradients in the flux points backup
  std::vector< std::vector< RealVector* > > m_cellGradFlxPntBackup;
  
  /// perturbation value
  RealVector m_eps;

  /// flags for each cell to tell whether its inner -divFD has been computed 
  std::vector< bool > m_cellFlags;
  
  /// unperturbed diffusive residuals for all cells
  std::vector< RealVector > m_unpertAllCellDiffRes;
  
  /// flux projection vectors in solution points for disc flux for a neighbor cell
  std::vector< std::vector< std::vector< RealVector > > > m_neighbCellFluxProjVects;
  
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

