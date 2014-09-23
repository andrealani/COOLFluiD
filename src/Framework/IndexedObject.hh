// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IndexedObject_hh
#define COOLFluiD_Framework_IndexedObject_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Framework/IndexList.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Object that can be indexed.
/// The object stores its index in a given domain to avoid searches in
/// connectivity tables.
/// It also stores a static data member to the global number list,
/// in order to ask it which ID should the object hold.
/// @author Tiago Quintino
/// @author Andrea Lani
/// @todo Why does this need to be a template?
template <class TYPE>
class IndexedObject {

  /// IndexList is made friend to have access to setLocalID
  friend class IndexList<TYPE>;

public:

  /// Default constructor without arguments
  IndexedObject();

  /// Copy constructor
  IndexedObject(const IndexedObject<TYPE>& idxObj);

  /// Default destructor
  ~IndexedObject();

  /// Assignment Operator
  IndexedObject<TYPE>& operator=(const IndexedObject<TYPE>& obj);

  /// An object can only be indexed by the appropriate IndexList
  /// @return if the object has been indexed.
  bool isIndexed() const
  {
    return (_localID != NO_INDEX && _globalID != NO_INDEX);
  }

  /// missing documentation
  /// @return the ID of this IndexedObject.
  CFuint getLocalID() const
  {
    cf_assert (hasLocalID());
    return _localID ;
  }

  /// missing documentation
  /// @return the ID of this IndexedObject.
  std::pair<CFuint, CFuint> getID() const
  {
    return std::pair<CFuint, CFuint>(getLocalID(),getGlobalID());
  }

  /// missing documentation
  ///@return the Global ID of this IndexedObject.
  CFuint getGlobalID() const
  {
    cf_assert (hasGlobalID());
    return _globalID ;
  }

  /// Set the Global ID of this IndexedObject.
  /// The local ID's are handed out by an IndexList,
  /// the Global ID's we can choose.
  void setGlobalID(const CFuint id)
  {
    cf_assert (id!=NO_INDEX);
    _globalID = id;
  }

  /// Check if a global index was set
  bool hasGlobalID () const
  {
    return _globalID != NO_INDEX;
  }

  /// check if a local index was set
  bool hasLocalID () const
  {
    return _localID != NO_INDEX;
  }

  /// Set the ID of this IndexedObject.
  /// @param id the ID to be the _localID
  void setLocalID(const CFuint id)
  {
    _localID = id;
  }

private:

  /// The index :-)
  CFuint _localID;

  /// The global index
  CFuint _globalID;

  /// This value means there was no index given to this object
  static const CFuint NO_INDEX;

}; // end of class IndexedObject

//////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
const CFuint IndexedObject<TYPE>::NO_INDEX = std::numeric_limits<CFuint>::max();

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "IndexedObject.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IndexedObject_hh
