// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NonCopyable
#define COOLFluiD_Common_NonCopyable

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Deriving from this class provides a clean and clear way to
/// show that a certain class is not copyable.
template < typename TYPE >
class NonCopyable {
public:

  /// Default inline constructor
  NonCopyable () {}

  /// Default inline destructor
  virtual ~NonCopyable () {}

private:

  /// private (non defined) copy constructor to prevent
  /// copying of the object
  NonCopyable (const NonCopyable & Source);

  /// private (non defined) assignment operator to prevent
  /// copy assignment of the object
  const NonCopyable & operator = (const NonCopyable & Source);

}; // end class NonCopyable

//////////////////////////////////////////////////////////////////////////////

    } // Common

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_NonCopyable
