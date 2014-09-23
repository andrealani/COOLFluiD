// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IndexList_hh
#define COOLFluiD_Framework_IndexList_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// The type of index's used by the IndexList
typedef CFuint IndexID;

/// This class represents IndexList of TYPE.
/// It defines a singleton for each TYPE. It is the implementation of the
/// Meyers Singleton pattern (Effective C++, Meyers 1998).
/// @author Tiago Quintino
/// @author Andrea Lani
template <class TYPE>
class IndexList {
public:

  /// the pointer to the TYPE
  typedef TYPE* TYPE_PTR;

  /// @return a reference to the sole instance of the Singleton IndexList.
  static IndexList<TYPE>& getList();

  /// @param ptr pointer to the object to index.
  /// @return a newly created IndexID for the object supplied.
  IndexID createID(TYPE_PTR ptr);

  /// @param ptr pointer to the object to index.
  /// @return a newly created IndexID for the object supplied.
  void removeID();

  /// Resets the IndexList.
  void reset()
  {
    _total = 0;
  }

  /// @return the size of the list
  CFuint size() const
  {
    return _total;
  }

private:

  /// Default destructor.
  ~IndexList();

  /// Default constructor without arguments
  IndexList();

  /// private copy constructor to insure singleton.
  IndexList(const IndexList&);

private:

  /// total number of pointers that are NOT CFNULL.
  CFuint _total;

}; // end of class IndexList

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "IndexList.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IndexList_hh
