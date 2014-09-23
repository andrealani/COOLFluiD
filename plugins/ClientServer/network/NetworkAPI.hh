#ifndef COOLFluiD_network_NetworkAPI_h
#define COOLFluiD_network_NetworkAPI_h

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro Network_API
/// @note build system defines network_EXPORTS when compiling Network files
#ifdef network_EXPORTS
#   define Network_API CF_EXPORT_API
#else
#   define Network_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_NetworkAPI_h
