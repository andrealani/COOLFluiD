#ifndef COOLFluiD_Physics_ArcJet_hh
#define COOLFluiD_Physics_ArcJet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a ArcJet flow physical model.
  namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ArcJet
 */
class ArcJetModule : public Environment::ModuleRegister<ArcJetModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ArcJet";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a ArcJet flow physical model.";
  }

}; // end ArcJetModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_ArcJet_hh
