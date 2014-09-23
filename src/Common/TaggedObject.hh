// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_TaggedObject_hh
#define COOLFluiD_Common_TaggedObject_hh

//////////////////////////////////////////////////////////////////////////////

#include <set>

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Object to which tags can be added
/// @author Tiago Quintino
class Common_API TaggedObject {

public: // type definitions

  typedef std::string Tag;

public: // functions

  /// Constructor with arguments
  TaggedObject();

  /// Default destructor
  ~TaggedObject();

  /// Checks if a tag is attached.
  /// @param tag tag to be checked for presence
  /// @return true if tag is attached
  bool hasTag(const Tag& tag) const;

  /// Checks if a tag is not attached.
  /// @param tag tag to be checked for absence
  /// @return true if tag is not attached
  bool hasNotTag(const Tag& tag) const;

  /// Attachs a tag to the object.
  /// If tag already exists it has no effect.
  /// @param tag tag to be attached
  void attachTag(const Tag& tag);

  /// Removes a tag from the object.
  /// If tag does not exist it has no effect.
  /// @param tag tag to be removed
  void removeTag(const Tag& tag);

  /// Prints the list of tags to the output stream separated with some given separator.
  /// @param out an output stream.
  /// @param sep the separator
  void print(std::ostream& out, const std::string& sep = std::string()) const;

private: // data

  /// storage of the tags
  std::set<Tag> m_tags;

}; // end of class TaggedObject

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_TaggedObject_hh
