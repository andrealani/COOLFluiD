#ifndef COOLFluiD_Numerics_FluctSplitNEQ_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    /// The classes that implement Fluctuation Splitting
    /// space method.
    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluctSplitNEQ
 */
class FluctSplitNEQModule : public Environment::ModuleRegister<FluctSplitNEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluctSplitNEQ";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implementes the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitNEQModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_hh
