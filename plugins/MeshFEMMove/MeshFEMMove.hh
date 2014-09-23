#ifndef COOLFluiD_Numerics_MeshFEMMove_hh
#define COOLFluiD_Numerics_MeshFEMMove_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

  /// The classes that implement a mesh adapter using a FEM method
  namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MeshFEMMove
 */
class MeshFEMMoveModule : public Environment::ModuleRegister<MeshFEMMoveModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MeshFEMMove";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a mesh adapter using a FEM method.";
  }

}; // end MeshFEMMoveModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_MeshFEMMove_hh

