#ifndef COOLFluiD_UFEM_UFEMQuadP1P1Cell_hh
#define COOLFluiD_UFEM_UFEMQuadP1P1Cell_hh

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro UFEM_API
/// @note build system defines UFEM_EXPORTS when compiling UFEM files
#ifdef UFEMQuadP1P1Cell_EXPORTS
#   define UFEMQuadP1P1Cell_API CF_EXPORT_API
#else
#   define UFEMQuadP1P1Cell_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace QuadP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

/// This class defines an Unsteady FEM solver module
class UFEMQuadP1P1Cell_API UFEMQuadP1P1CellPlugin : public Environment::ModuleRegister< UFEMQuadP1P1CellPlugin > {
public:

  /// Static function that returns the module name.
  /// @return name of the module
  static std::string getModuleName() { return "UFEMQuadP1P1Cell"; }

  /// Static function that returns the description of the module.
  /// @return descripton of the module
  static std::string getModuleDescription() { return "This module is a Unsteady FEM solver." ; }

}; // end UFEMPlugin

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace QuadP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_UFEMQuadP1P1Cell_hh

