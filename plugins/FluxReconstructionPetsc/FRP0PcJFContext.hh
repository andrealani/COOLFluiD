// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionPetsc_FRP0PcJFContext_hh
#define COOLFluiD_FluxReconstructionPetsc_FRP0PcJFContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // for Mat, Vec, KSP types
#include "Framework/DataStorage.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    class JFContext;

//////////////////////////////////////////////////////////////////////////////

/**
 * Context data for the two-level p-multigrid FR shell preconditioner.
 *
 * Supports two block modes:
 *   - ElementBlock: full (nSolPts*nEqs)^2 element-diagonal blocks (smoother)
 *                   + Galerkin-projected nEqs^2 P0 blocks (coarse)
 *   - PointBlock:   per-DOF nEqs^2 diagonal sub-blocks (smoother)
 *                   + averaged nEqs^2 P0 blocks (coarse)
 *
 * P0 blocks are derived from the assembled Jacobian via Galerkin projection
 * (no separate P0 matrix, no inner KSP).
 *
 * @author Rayan Dhib
 */
class FRP0PcJFContext {
public:

  /// Constructor
  FRP0PcJFContext() : pJFC(CFNULL), nEqs(0), nDim(0), nUpdatableCells(0),
    nSolPtsPerCell(0), useElementBlocks(true), smootherOmega(1.0),
    useP0ILU(false), p0Mat(CFNULL), p0RhsVec(CFNULL), p0SolVec(CFNULL),
    p0Ksp(CFNULL) {}

  /// For each updatable cell: list of indices into the PETSc updatable state vector
  /// cellToUpStateIdx[iCell][iSol] = updatable state index for sol pt iSol in cell iCell
  std::vector<std::vector<CFuint> > cellToUpStateIdx;

  /// Number of equations
  CFuint nEqs;

  /// Spatial dimension
  CFuint nDim;

  /// Number of updatable cells
  CFuint nUpdatableCells;

  /// Pointer to JFContext
  JFContext* pJFC;

  /// Number of sol pts per cell (uniform across all cells)
  CFuint nSolPtsPerCell;

  // ---- Mode flag ----

  /// true = ElementBlock mode, false = PointBlock mode
  bool useElementBlocks;

  // ---- Smoother data ----

  /// ElementBlock mode: per-cell inverted (nSolPts*nEqs x nSolPts*nEqs) blocks
  std::vector<RealMatrix> invElemBlocks;

  /// PointBlock mode: per-DOF inverted (nEqs x nEqs) diagonal sub-blocks
  std::vector<RealMatrix> smootherInvBlocks;

  /// Relaxation parameter for the smoother (default 1.0)
  CFreal smootherOmega;

  // ---- P0 coarse correction data (both modes) ----

  /// Per-cell inverted (nEqs x nEqs) P0 blocks from Galerkin projection
  /// Used in BlockDiag mode (useP0ILU == false)
  std::vector<RealMatrix> p0InvBlocks;

  /// P0 restricted residual buffer, size nCells*nEqs
  std::vector<CFreal> p0Residual;

  /// P0 correction buffer, size nCells*nEqs
  std::vector<CFreal> p0Correction;

  // ---- P0 ILU coarse solve data (when useP0ILU == true) ----

  /// true = use face-coupled P0 sparse matrix with inner ILU-GMRES solve
  /// false = use element-diagonal P0 block inversion (legacy)
  bool useP0ILU;

  /// Cell-to-cell adjacency: cellNeighbors[iUpdCell] = updatable cell indices of face neighbors
  std::vector<std::vector<CFuint> > cellNeighbors;

  /// Mapping from TRS cell index to updatable cell index (-1 if not updatable)
  std::vector<CFint> cellToUpdIdx;

  /// PETSc SeqBAIJ matrix for the face-coupled P0 coarse level
  Mat p0Mat;

  /// PETSc Vecs for P0 coarse level rhs and solution
  Vec p0RhsVec;
  Vec p0SolVec;

  /// Inner KSP for P0 coarse solve (GMRES + ILU)
  KSP p0Ksp;

}; // end of class FRP0PcJFContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_FRP0PcJFContext_hh
