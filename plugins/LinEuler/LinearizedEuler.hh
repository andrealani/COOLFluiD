#ifndef COOLFluiD_Physics_LinearizedEuler_hh
#define COOLFluiD_Physics_LinearizedEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Linearized Euler
 * @author Lilla Edit Koloszar
 * @author Nadege Villedieu
 * @author Tomas Kopacek
 */
class LinearizedEulerModule : public Environment::ModuleRegister<LinearizedEulerModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "LinearizedEuler";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Linearized Euler flow physical model.";
  }

}; // end LinearizedEulerModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_LinEuler_hh
