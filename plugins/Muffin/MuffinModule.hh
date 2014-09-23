#ifndef COOLFluiD_Muffin_MuffinModule_hh
#define COOLFluiD_Muffin_MuffinModule_hh

#include "Environment/ModuleRegister.hh"

namespace COOLFluiD {

  /// The classes that implement RDS space discretization method.
  /// Solves incompressible flow with temperature and MITReM equations.
  namespace Muffin {

/// This class defines the Module Muffin
class MuffinModule : public Environment::ModuleRegister< MuffinModule > {
public:

  /// Static function that returns the module name. Must be implemented for the
  /// ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {
    return "Muffin";
  }

  /// Static function that returns the description of the module. Must be
  /// implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() {
    return "This module interfaces the Muffin solver.";
  }

}; // end MuffinModule


  }  // end namespace Muffin
}  // end namespace COOLFluiD

#endif // COOLFluiD_Muffin_MuffinModule_hh

