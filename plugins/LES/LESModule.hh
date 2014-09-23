#ifndef COOLFluiD_LES_LESModule_hh
#define COOLFluiD_LES_LESModule_hh

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace LES {

//////////////////////////////////////////////////////////////////////////////

/// Class defining the module LESModule
class LESModule : public Environment::ModuleRegister< LESModule > {
public:

  /// Static function that returns the module name (must be implemented for
  /// the ModuleRegister template)
  static std::string getModuleName()
  {
    return "LESModule";
  }

  /// Static function that returns the description of the module (must be
  /// implemented for the ModuleRegister template)
  static std::string getModuleDescription()
  {
    return "This module interfaces the LES particle tracking library";
  }

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace LES
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFLUID_LES_LESModule_hh

