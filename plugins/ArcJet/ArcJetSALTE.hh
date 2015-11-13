#ifndef COOLFluiD_Physics_ArcJetSALTE_hh
#define COOLFluiD_Physics_ArcJetSALTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a ArcJetSALTE flow physical model.
  namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ArcJetSALTE
 */
class ArcJetSALTEModule : public Environment::ModuleRegister<ArcJetSALTEModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ArcJetSALTE";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a ArcJetSALTE flow physical model.";
  }

}; // end ArcJetSALTEModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_ArcJetSALTE_hh
