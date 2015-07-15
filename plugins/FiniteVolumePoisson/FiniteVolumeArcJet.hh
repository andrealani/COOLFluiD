#ifndef COOLFluiD_Numerics_FiniteVolumeArcJet_hh
#define COOLFluiD_Numerics_FiniteVolumeArcJet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /// Classes concerning Inductively Coupled Plasma (FiniteVolumeArcJet)
    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeArcJet
 */
class FiniteVolumeArcJetModule : public Environment::ModuleRegister<FiniteVolumeArcJetModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeArcJet";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module concerns Inductively Coupled Plasma : FiniteVolumeArcJet.";
  }

}; // end FiniteVolumeArcJetModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeArcJet_hh
