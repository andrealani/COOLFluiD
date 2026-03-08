// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionPetsc_FRBlockJacobiPcJFContext_hh
#define COOLFluiD_FluxReconstructionPetsc_FRBlockJacobiPcJFContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataStorage.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    class JFContext;

//////////////////////////////////////////////////////////////////////////////

/**
 * Context data for the FR element block-Jacobi shell preconditioner.
 *
 * Stores per-cell inverted diagonal blocks and the cell-to-state mapping
 * needed by the PCShellApply callback.
 *
 * @author Rayan Dhib
 */
class FRBlockJacobiPcJFContext {
public: // functions

  /// Constructor
  FRBlockJacobiPcJFContext() : pJFC(CFNULL), nEqs(0), nUpdatableCells(0) {}

  /// For each updatable cell: list of updatable state indices (position in PETSc vector / nEqs)
  /// cellToUpStateIdx[iCell][iSol] = sequential updatable index for sol pt iSol in cell iCell
  std::vector<std::vector<CFuint> > cellToUpStateIdx;

  /// Per-cell inverted element diagonal blocks
  /// invElemBlocks[iCell] has size (nSolPts*nEqs x nSolPts*nEqs)
  std::vector<RealMatrix> invElemBlocks;

  /// Number of equations
  CFuint nEqs;

  /// Number of updatable cells
  CFuint nUpdatableCells;

  /// Pointer to JFContext (provides states, rhs, updateCoeff, spaceMethod, petscData)
  JFContext* pJFC;

}; // end of class FRBlockJacobiPcJFContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_FRBlockJacobiPcJFContext_hh
