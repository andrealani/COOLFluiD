#ifndef COOLFluiD_PLaS_PLaSModule_hh
#define COOLFluiD_PLaS_PLaSModule_hh

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// Class defining the module PLaSModule
class PLaSModule : public Environment::ModuleRegister< PLaSModule > {
public:

  /// Static function that returns the module name (must be implemented for
  /// the ModuleRegister template)
  static std::string getModuleName()
  {
    return "PLaSModule";
  }

  /// Static function that returns the description of the module (must be
  /// implemented for the ModuleRegister template)
  static std::string getModuleDescription()
  {
    return "This module interfaces the PLaS particle tracking library";
  }

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFLUID_PLaS_PLaSModule_hh

