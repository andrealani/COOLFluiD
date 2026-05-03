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
 * Additive mode:
 *   Y = omega * B^{-1} * X  +  P * D0^{-1} * R * X
 *
 * Multiplicative mode (default):
 *   y = omega * B^{-1} * X                (smoother)
 *   d = X - A * y                          (defect)
 *   Y = y + P * D0^{-1} * R * d           (coarse correction on defect)
 *
 * Block modes (config BlockMode):
 *   - ElementBlock: full (nSolPts*nEqs)^2 element-diagonal blocks
 *   - PointBlock:   per-DOF (nEqs)^2 diagonal sub-blocks
 *
 * Coarse solve modes (config CoarseSolveType):
 *   - FaceCoupled: Galerkin-projected face-coupled P0 sparse matrix, direct LU
 *   - BlockDiag: element-local nEqs x nEqs block inversion (no face coupling)
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

  /// P0 coarse solve type: "FaceCoupled" or "BlockDiag"
  std::string _coarseSolveType;

  /// Use multiplicative (defect-based) coarse correction (default false)
  bool _multiplicative;

}; // end of class FRP0Preconditioner

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_FRP0Preconditioner_hh
