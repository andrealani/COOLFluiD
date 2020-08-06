#ifndef COOLFluiD_Physics_GammaAlpha_hh
#define COOLFluiD_Physics_GammaAlpha_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a GammaAlpha flow physical model.
  namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module GammaAlpha
 *
 * @author Ray Vandenhoeck
  */


class GammaAlphaModule : public Environment::ModuleRegister<GammaAlphaModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "GammaAlpha";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a GammaAlpha flow physical model.";
  }

}; // end GammaAlphaModule

//////////////////////////////////////////////////////////////////////////////

    } // end namespace GammaAlpha

  } // end namespace Physics

} // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_GammaAlpha_hh
