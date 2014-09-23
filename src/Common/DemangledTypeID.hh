// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_DemangledTypeID_hh
#define COOLFluiD_Common_DemangledTypeID_hh

/// @note This header should be included by including COOLFluiD.hh instead.

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Function to demangle the return of typeid()
Common_API std::string demangle (const char* type);

#define DEMANGLED_TYPEID(a) COOLFluiD::Common::demangle(typeid(a).name())

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_DemangledTypeID_hh
