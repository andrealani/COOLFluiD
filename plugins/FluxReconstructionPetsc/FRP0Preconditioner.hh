// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionPetsc_FRP0Preconditioner_hh
#define COOLFluiD_FluxReconstructionPetsc_FRP0Preconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/ShellPreconditioner.hh"
#include "FluxReconstructionPetsc/FRP0PcJFContext.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { class MatrixInverter; }

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * Two-level additive p-multigrid shell preconditioner for JFNK in FR solvers.
 *
 * Apply formula:
 *   Y = omega * B^{-1} * X  +  P * D0^{-1} * R * X
 *       (smoother)            (P0 coarse correction)
 *
 * Two block modes (config option BlockMode):
 *   - ElementBlock: smoother uses full (nSolPts*nEqs)^2 element-diagonal blocks
 *   - PointBlock:   smoother uses per-DOF (nEqs)^2 diagonal sub-blocks
 *
 * P0 blocks are derived from the assembled Jacobian via Galerkin projection
 * (R * B * P). Two coarse solve modes (config option CoarseSolveType):
 *   - BlockDiag (default): element-local nEqs x nEqs block inversion
 *   - ILU: face-coupled P0 sparse matrix with inner GMRES + ILU(0) solve
 *
 * @author Rayan Dhib
 */
class FRP0Preconditioner : public ShellPreconditioner {
public:

  /// Define config options
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  FRP0Preconditioner(const std::string& name);

  /// Destructor
  ~FRP0Preconditioner();

  /// Set up: build cell-to-state mapping, allocate blocks, register PCShellApply
  virtual void setPreconditioner();

  /// Assemble Jacobian, extract and invert smoother + P0 blocks
  virtual void computeBeforeSolving();

  /// Cleanup: zero all block storage
  virtual void computeAfterSolving();

  /// Returns needed data sockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for updateCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for volumes (optional)
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// P0 preconditioner context (passed to PCShellApply)
  FRP0PcJFContext _pcc;

  /// Matrix inverter for smoother blocks (ElementBlock: nSolPts*nEqs, PointBlock: nEqs)
  std::auto_ptr<MathTools::MatrixInverter> _inverter;

  /// Matrix inverter for P0 blocks (always nEqs x nEqs)
  std::auto_ptr<MathTools::MatrixInverter> _p0Inverter;

  /// Smoother relaxation parameter (configurable via SmootherOmega)
  CFreal _smootherOmega;

  /// Block mode: "ElementBlock" or "PointBlock"
  std::string _blockMode;

  /// Flag: compute element blocks directly without PETSc matrix
  bool _directBlocks;

  /// P0 coarse solve type: "ILU" or "BlockDiag"
  std::string _coarseSolveType;

  /// Inner KSP relative tolerance for P0 coarse solve (default 0.1)
  CFreal _coarseKSPTol;

  /// Inner KSP max iterations for P0 coarse solve (default 5)
  CFuint _coarseMaxIter;

  /// Direct element-diagonal block storage (indexed by cell TRS-local index)
  std::vector<RealMatrix> m_cellBlocks;

}; // end of class FRP0Preconditioner

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_FRP0Preconditioner_hh
