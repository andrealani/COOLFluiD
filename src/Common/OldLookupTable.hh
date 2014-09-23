// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_LookUpTable_hh
#define COOLFluiD_Common_LookUpTable_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents table where values can be looked up
/// and if not found, are interpolated with the nearest ones.
/// @todo Add trait class for interpolation algorithm
/// @author Tiago Quintino
template <class KEY, class VALUE>
class LookUpTable : public CFMap<KEY,VALUE> {
public:

  typedef KEY   key_type;
  typedef VALUE value_type;
  typedef LookUpTable<key_type,value_type> lookuptable_type;

public:

  /// Constructor
  LookUpTable(size_t maxSize = 0);

  /// Default destructor
  ~LookUpTable();

  /// Gets the interpolated VALUE of the supplied KEY.
  /// @param key entry KEY to be interpolated
  VALUE get(const KEY& key);

  /// Outputs the table to the supplied stream
  /// @param outStream the stream to output the table
  void saveToStream(std::ostream& outStream);

  /// Inputs the table from the supplied stream
  /// @param inStream the stream to input the table from.
  void readFromStream(std::istream& inStream);

}; // end of class LookUpTable

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "OldLookupTable.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_LookUpTable_ci
