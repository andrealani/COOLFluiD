#ifndef COOLFluiD_Numerics_MatrixStabilityMethodWriter_hh
#define COOLFluiD_Numerics_MatrixStabilityMethodWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

  /// The classes that implement a tool that writes the matrix corresponding.to a spatial discretization
  /// to a file. This matrix can then be used for a stability analysis based on the matrix method.
  /// @see "Numerical computation of internal and external flows, Volume 1, C. Hirsch, Chapter 10
  /// @note should be used in combination with a Linear Advection/Linear Diffusion physical model
  namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MatrixStabilityMethodWriter
 */
class MatrixStabilityMethodWriterModule : public Environment::ModuleRegister<MatrixStabilityMethodWriterModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MatrixStabilityMethodWriter";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a tool that writes the matrix corresponding to the spatial method for a certain test case. To be used for matrix method for stability analysis";
  }

}; // end MatrixStabilityMethodWriterModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_MatrixStabilityMethodWriter_hh
