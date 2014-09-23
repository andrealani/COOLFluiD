#ifndef COOLFluiD_UFEM_UFEM_hh
#define COOLFluiD_UFEM_UFEM_hh

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro UFEM_API
/// @note build system defines UFEM_EXPORTS when compiling UFEM files
#ifdef UFEM_EXPORTS
#   define UFEM_API CF_EXPORT_API
#   define UFEM_TEMPLATE
#else
#   define UFEM_API CF_IMPORT_API
#   define UFEM_TEMPLATE CF_TEMPLATE_EXTERN
#endif

/// Macro defining the factories in Environment
#ifdef  CF_HAVE_CXX_EXPLICIT_TEMPLATES
#   define UFEM_Factory(__fac__) namespace COOLFluiD { namespace Environment { CF_TEMPLATE_EXTERN template class UFEM_API COOLFluiD::Environment::Factory<COOLFluiD::UFEM::__fac__>; } }
#else
#   define UFEM_Factory(__fac__)
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class defines an Unstructured FEM solver module
class UFEM_API UFEMPlugin : public Environment::ModuleRegister< UFEMPlugin > {
public:

  /// Static function that returns the module name.
  /// @return name of the module
  static std::string getModuleName() { return "UFEM"; }

  /// Static function that returns the description of the module.
  /// @return descripton of the module
  static std::string getModuleDescription() { return "This module is a Unsteady FEM solver." ; }

}; // end UFEMPlugin

//////////////////////////////////////////////////////////////////////////////

  }  // namespace UFEM
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_UFEM_hh

