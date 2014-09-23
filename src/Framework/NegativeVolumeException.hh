// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NegativeVolumeException_hh
#define COOLFluiD_Framework_NegativeVolumeException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the exception thrown if a certain path
/// to a directory does not exist
/// @author Tiago Quintino
class Framework_API NegativeVolumeException : public Common::Exception {
public:

  /// Constructor
  NegativeVolumeException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "NegativeVolumeException") {}

  /// Copy constructor
  NegativeVolumeException ( const NegativeVolumeException& e) throw () : Exception(e) {}

}; // end of class NegativeVolumeException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NegativeVolumeException_hh
