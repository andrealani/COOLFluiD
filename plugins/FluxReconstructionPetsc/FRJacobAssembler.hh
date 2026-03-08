// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionPetsc_FRJacobAssembler_hh
#define COOLFluiD_FluxReconstructionPetsc_FRJacobAssembler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/ShellPreconditioner.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * Assembly-only shell preconditioner for JFNK in FR solvers.
 *
 * This preconditioner assembles the full Jacobian into the PETSc
 * preconditioner matrix (precondMat) and does nothing else.
 * PETSc's built-in PC (e.g. ILU, ASM, BJACOBI) operates directly
 * on precondMat.
 *
 * Unlike FRBlockJacobi or FRP0Precond, this class does NOT call
 * PCShellSetApply(), so the PCType set in the config (PCILU, etc.)
 * is preserved and used by PETSc natively.
 *
 * Config: DifferentPreconditionerMatrix = true
 *         UseBlockPreconditioner = false  (ParBAIJ path)
 *         PCType = PCILU (or PCASM, PCBJACOBI, etc.)
 *         ShellPreconditioner = FRJacobAssembler
 *
 * @author Rayan Dhib
 */
class FRJacobAssembler : public ShellPreconditioner {
public:

  /// Constructor
  FRJacobAssembler(const std::string& name);

  /// Destructor
  ~FRJacobAssembler();

  /// No-op: we don't register a PCShellApply callback,
  /// so PETSc's built-in PC (ILU, ASM, etc.) is used instead.
  virtual void setPreconditioner();

  /// Assemble the full Jacobian into the preconditioner matrix.
  virtual void computeBeforeSolving();

  /// No-op
  virtual void computeAfterSolving();

  /// Returns needed data sockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for updateCoeff (needed to backup/restore during assembly)
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

}; // end of class FRJacobAssembler

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_FRJacobAssembler_hh
