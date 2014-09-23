#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement MeshAdapterSpringAnalogy mesh adapter.
  namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MeshAdapterSpringAnalogy
 */
class MeshAdapterSpringAnalogyModule : public Environment::ModuleRegister<MeshAdapterSpringAnalogyModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MeshAdapterSpringAnalogy";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a mesh adapter based on spring analogy.";
  }

}; // end MeshAdapterSpringAnalogyModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_MeshAdapterSpringAnalogy_hh
