#ifndef COOLFluiD_Numerics_FiniteVolumePoissonNEQ_hh
#define COOLFluiD_Numerics_FiniteVolumePoissonNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /// Classes concerning PoissonNEQ (FiniteVolumePoissonNEQ)
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumePoissonNEQ
 */
class FiniteVolumePoissonNEQModule : public Environment::ModuleRegister<FiniteVolumePoissonNEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumePoissonNEQ";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module concerns PoissonNEQ Equation solver : FiniteVolumePoissonNEQ.";
  }

}; // end FiniteVolumePoissonNEQModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumePoissonNEQ_hh
