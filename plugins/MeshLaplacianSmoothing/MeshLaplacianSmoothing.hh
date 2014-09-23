#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

  /// The classes that implement a mesh adapter using a Laplacian smoothing
  namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MeshLaplacianSmoothing
 */
class MeshLaplacianSmoothingModule : public Environment::ModuleRegister<MeshLaplacianSmoothingModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MeshLaplacianSmoothing";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements  a mesh adapter using a Laplacian smoothing.";
  }

}; // end MeshLaplacianSmoothingModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_MeshLaplacianSmoothing_hh

