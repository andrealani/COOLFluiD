#ifndef COOLFluiD_Numerics_RKRD_hh
#define COOLFluiD_Numerics_RKRD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement a general RKRD time stepper.
  namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module RKRD
class RKRDModule : public Environment::ModuleRegister<RKRDModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "RKRD";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements a general RKRD time stepper.";
  }

}; // end RKRDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_RKRD_hh
