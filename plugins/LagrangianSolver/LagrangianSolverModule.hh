#ifndef COOLFluiD_LagrangianSolver_LagrangianSolverModule_hh
#define COOLFluiD_LagrangianSolver_LagrangianSolverModule_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module RadiativeTransfer
 */
class LagrangianSolverModule : public Environment::ModuleRegister<LagrangianSolverModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "RadiativeTransfer";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the Radiative Heat Transfer.";
  }

}; // end RadiativeTransfer

//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_RadiativeTransfer_RadiativeTransfer_hh
