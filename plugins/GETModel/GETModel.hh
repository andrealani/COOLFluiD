#ifndef COOLFluiD_Numerics_GETModel_hh
#define COOLFluiD_Numerics_GETModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement an DomainModel which provides
  /// an CAD description of the computational domain
  namespace GETModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module GETModel
 *
 * @author Tiago Quintino
 */
class GETModelModule : public Environment::ModuleRegister<GETModelModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "GETModel";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an analytical description of the computational domain";
  }

}; // end GETModelModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace GETModel

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_GETModel_hh

