#ifndef COOLFluiD_DiscontGalerkin_hh
#define COOLFluiD_DiscontGalerkin_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement computation by discontinuous Galerkin method
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module DiscontGalerkin
 */
class DiscontGalerkinModule : public Environment::ModuleRegister<DiscontGalerkinModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "DiscontGalerkin";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements  computation by discontinuous Galerkin method.";
  }
}; // end DiscontGalerkinModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace DiscontGalerkin

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_DiscontGalerkin_hh

