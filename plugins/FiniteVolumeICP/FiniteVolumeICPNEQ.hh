#ifndef COOLFluiD_Numerics_FiniteVolumeICPNEQ_hh
#define COOLFluiD_Numerics_FiniteVolumeICPNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /// Classes concerning Inductively Coupled Plasma (FiniteVolumeICPNEQ)
    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeICPNEQ
 */
class FiniteVolumeICPNEQModule : public Environment::ModuleRegister<FiniteVolumeICPNEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeICPNEQ";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module concerns Inductively Coupled Plasma : FiniteVolumeICPNEQ.";
  }

}; // end FiniteVolumeICPNEQModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeICPNEQ_hh
