// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_StringOps_hh
#define COOLFluiD_StringOps_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonInstantiable.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Operations on std::string
/// @author Tiago Quintino
class Common_API StringOps : public Common::NonInstantiable<StringOps> {

public:

  /// Joins the array of parts, using the delimiter.
  /// @param delim Missing documentation
  /// @param parts Missing Documentation
  /// @param nparts Missing documentation
  /// @param out string to modify
  static void join (const std::string& delim, std::string parts[], unsigned int nparts, std::string& out);

  /// Modifies string to lowercase.
  /// @param out string to modify
  static void toLower (std::string& out);

  /// Modifies string to uppercase.
  /// @param out string to modify
  static void toUpper (std::string& out);

  /// Substitute the rhs for the lhs wherever it appears in the string
  /// @param lhs partial std::string to match
  /// @param rhs partial std::string to substitute
  /// @param out string to modify
  static void subst (const std::string& lhs, const std::string& rhs, std::string& out);

  /// Removes the leading whitespace.
  /// @param out string to modify
  static void trimFront (std::string& out);

  /// Removes the trailing whitespaces.
  /// @param out string to modify
  static void trimRear (std::string& out);

  /// Removes the trailing and leading whitespace.
  /// @param out string to modify
  static void trim (std::string& out);

  /// Removes the trailing and leading whitespace.
  /// @param out string to modify
  static void trim2 (std::string& out);
  
  /// Splits the given std::string on whitespace, and puts the words into the given vector.
  /// @param in string to split
  /// @return a vector with the separated std::string's
  static std::vector<std::string> getWords (const std::string& in);

  /// Splits the given std::string on the passed characters, and puts the words into the given vector.
  /// @param in string to split
  /// @param sep the separator character
  /// @return a vector with the separated std::string's
  static std::vector<std::string> getWords  (const std::string& in, const CFchar sep );

  /// Returns whether this std::string starts with the given std::string.
  /// @param in string to split
  /// @param str missing documentation
  /// @return missing documentation
  static bool startsWith (const std::string& in, const std::string& str);
  
  /// Returns whether this std::string ends with the given std::string.
  /// @param in string to split
  /// @param str missing documentation
  /// @return missing documentation
  static bool endsWith (const std::string& in, const std::string& str);

  /// Converts to std::string
  /// Don't use this to convert to a CFchar, use c_str for that.
  /// Typical use is to convert to numbers.
  /// @param str string to convert from
  /// @return converter type
  template <class T>
  static std::string to_str (T v)
  {
    std::ostringstream oss;
    oss << v;
    return oss.str();
  }

  /// Converts from std::string
  /// Don't use this to convert to a CFchar, use c_str for that.
  /// Typical use is to convert to numbers.
  /// @param str string to convert from
  /// @return converter type
  template <class T>
  static T from_str (const std::string& str)
  {
    T v;
    if (str.length() > 0)
    {
      std::istringstream iss(str.c_str());
      iss >> v;
    }
    else
    {
      // pretty much everything has an empty constuctor
      v = T();
    }
    return v;
  }

}; // class StringOps

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_StringOps_hh
