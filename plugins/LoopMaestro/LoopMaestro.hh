#ifndef COOLFluiD_Numerics_LoopMaestro_hh
#define COOLFluiD_Numerics_LoopMaestro_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement a Maestro that loops on subsystems
  namespace LoopMaestro {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module LoopMaestro
class LoopMaestroPlugin : public Environment::ModuleRegister<LoopMaestroPlugin> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() { return "LoopMaestro"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements a Maestro that loops on subsystems.";
  }

}; // end LoopMaestroPlugin

//////////////////////////////////////////////////////////////////////////////

} // namespace LoopMaestro
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_LoopMaestro_hh

