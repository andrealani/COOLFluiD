#ifndef COOLFluiD_FluxReconstructionEntropyMHD_hh
#define COOLFluiD_FluxReconstructionEntropyMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionEntropyMHD
 */
class FluxReconstructionEntropyMHDModule : public Environment::ModuleRegister<FluxReconstructionEntropyMHDModule> {
public:

  static std::string getModuleName()
  {
    return "FluxReconstructionEntropyMHD";
  }

  static std::string getModuleDescription()
  {
    return "This module implements FR-specific BCs and source terms for the EntropyMHD model.";
  }

}; // end FluxReconstructionEntropyMHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionEntropyMHD_hh
