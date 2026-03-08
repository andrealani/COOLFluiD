// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionPetsc_FRBlockJacobiPreconditioner_hh
#define COOLFluiD_FluxReconstructionPetsc_FRBlockJacobiPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/ShellPreconditioner.hh"
#include "FluxReconstructionPetsc/FRBlockJacobiPcJFContext.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { class MatrixInverter; }

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * FR element block-Jacobi shell preconditioner for JFNK solvers.
 *
 * This preconditioner assembles the full Jacobian into the PETSc
 * preconditioner matrix, extracts per-cell diagonal blocks
 * (nSolPts*nEqs x nSolPts*nEqs), inverts them, and applies the
 * inverse as a block-Jacobi preconditioner during GMRES iterations.
 *
 * Unlike the FRP0Preconditioner PointBlock mode, this captures:
 * - Inter-equation coupling (dF_rho/d(rhoU), etc.)
 * - Intra-cell coupling between solution points
 * - Per-sol-pt time diagonal variation
 *
 * Memory: N_cells * (nSolPts*nEqs)^2 * 8 bytes per cell.
 * For 2D P3 quad, 4 eqs: 2328 cells * 64^2 * 8 = ~76 MB.
 *
 * @author Rayan Dhib
 */
class FRBlockJacobiPreconditioner : public ShellPreconditioner {
public:

  /// Constructor
  FRBlockJacobiPreconditioner(const std::string& name);

  /// Destructor
  ~FRBlockJacobiPreconditioner();

  /// Set up the preconditioner: build cell-to-state mapping, register PCShellApply
  virtual void setPreconditioner();

  /// Assemble Jacobian, extract per-cell diagonal blocks, invert them
  virtual void computeBeforeSolving();

  /// Cleanup: zero out inverted blocks
  virtual void computeAfterSolving();

  /// Returns needed data sockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for updateCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// FR block-Jacobi preconditioner context
  FRBlockJacobiPcJFContext _pcc;

  /// Matrix inverter (LU-based, allocated in setPreconditioner for the element block size)
  std::auto_ptr<MathTools::MatrixInverter> _inverter;

}; // end of class FRBlockJacobiPreconditioner

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_FRBlockJacobiPreconditioner_hh
