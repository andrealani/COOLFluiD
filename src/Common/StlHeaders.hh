// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_StlHeaders_hh
#define COOLFluiD_Common_StlHeaders_hh

/// @note This header should be included by including COOLFluiD.hh instead.

//////////////////////////////////////////////////////////////////////////////

/// @note This file must not include iostream because it slows compilation.
///       Files that need iostream should include it directly.

#include <memory>    // for std::auto_ptr
#include <string>    // for std::string
#include <sstream>   // string stream
#include <vector>    // for std::vector
#include <map>       // for std::map
#include <cmath>     // all sorts of mathematical functions and constants
#include <algorithm> // all sorts of algorithms on stl containers
#include <utility>   // for stl::pair
#include <typeinfo>  // for typeid
#include <complex>   // for complex numbers and complex functions

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_StlHeaders_hh
