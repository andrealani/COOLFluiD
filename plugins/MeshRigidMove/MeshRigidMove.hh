#ifndef COOLFluiD_Numerics_MeshRigidMove_hh
#define COOLFluiD_Numerics_MeshRigidMove_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement MeshRigidMove mesh adapter.
  namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MeshRigidMove
 */
class MeshRigidMoveModule : public Environment::ModuleRegister<MeshRigidMoveModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MeshRigidMove";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a rigid move mesh adapter.";
  }

}; // end MeshRigidMoveModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_MeshRigidMove_hh
