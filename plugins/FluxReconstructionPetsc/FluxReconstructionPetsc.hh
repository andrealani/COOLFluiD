#ifndef COOLFluiD_FluxReconstructionPetsc_hh
#define COOLFluiD_FluxReconstructionPetsc_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionPetsc
 *
 * Bridge plugin providing FR-specific PETSc preconditioners and Jacobian
 * assemblers for JFNK solvers. Depends on both PetscI and
 * FluxReconstructionMethod.
 *
 * Note: The preconditioner/assembler classes remain in the COOLFluiD::Petsc
 * namespace because they inherit from Petsc::ShellPreconditioner and use
 * Petsc types extensively. Only the module class lives in
 * COOLFluiD::FluxReconstructionMethod, following the FluxReconstructionMHD
 * pattern.
 */
class FluxReconstructionPetscModule :
    public Environment::ModuleRegister<FluxReconstructionPetscModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionPetsc";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return description of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements FR-specific PETSc preconditioners and "
           "Jacobian assemblers for JFNK solvers.";
  }

}; // end FluxReconstructionPetscModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionPetsc_hh
