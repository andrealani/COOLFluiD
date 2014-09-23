#ifndef COOLFluiD_IO_TecplotWriterNavierStokes_hh
#define COOLFluiD_IO_TecplotWriterNavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Tecplot writer
  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module TecplotWriterNavierStokes
 */
class TecplotWriterNavierStokesModule : public Environment::ModuleRegister<TecplotWriterNavierStokesModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "TecplotWriterNavierStokes";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the Tecplot writer specific for Navier-Stokes";
  }

}; // end TecplotWriterModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_TecplotWriterNavierStokes_hh
