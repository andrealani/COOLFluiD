// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NamedObject_hh
#define COOLFluiD_Common_NamedObject_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Object that can be named.
/// @author Tiago Quintino
class Common_API NamedObject {

public:

  /// Constructor with arguments
  explicit NamedObject(const std::string& name = std::string());
  /// Default destructor
  virtual ~NamedObject();

  /// Gets the name of the object.
  /// @return std::string with the object name.
  std::string getName() const { return m_name; }

protected: // functions

  /// Sets the object name
  /// @param name std::string with object name
  void setName(const std::string& name) {  m_name = name;  }

private: // data

  /// The object name stored as a std::string
  std::string m_name;

}; // end of class NamedObject

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_NamedObject_hh
