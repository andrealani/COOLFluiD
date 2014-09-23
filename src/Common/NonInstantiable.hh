// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NonInstantiable_hh
#define COOLFluiD_Common_NonInstantiable_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Derive from this class if you want a class that is not instantiable.
template < class TYPE >
class NonInstantiable {
private:
    NonInstantiable ();
};

//////////////////////////////////////////////////////////////////////////////

    } // namespace Common
} // namespace COOLFluiD

#endif // COOLFluiD_Common_NonInstantiable_hh
