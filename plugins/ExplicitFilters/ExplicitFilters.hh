#ifndef COOLFluiD_Numerics_ExplicitFilters_hh
#define COOLFluiD_Numerics_ExplicitFilters_hh

//////////////////////////////////////////////////////////////////////////////

/*
 * @todo  Check which headers are required everywhere
 * @todo  Implement way to gather stencil  
 * 					- make geobuilder
 * 					- look at ComputeStencil in FiniteVolume for example
 */

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /** 
     * The classes that implement the ExplicitFilters method
     */
    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ExplicitFilters
 * 
 * @author Willem Deconinck
 */
class ExplicitFiltersModule : public Environment::ModuleRegister<ExplicitFiltersModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ExplicitFilters";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the ExplicitFilters method.";
  }

}; // end ExplicitFiltersModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_ExplicitFilters_hh

