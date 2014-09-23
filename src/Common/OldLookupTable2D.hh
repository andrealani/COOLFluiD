// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_LookUpTable2D_hh
#define COOLFluiD_Common_LookUpTable2D_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "LookUpTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents table where values can be looked up
/// and if not found, are interpolated with the nearest ones.
/// @todo Add trait class for interpolation algorithm
/// @author Tiago Quintino
template <class KEY1, class KEY2, class VALUE>
class LookUpTable2D {

  typedef KEY1  key1_type;
  typedef KEY2  key2_type;
  typedef VALUE value_type;
  typedef LookUpTable<KEY2,VALUE> intable_type;
  typedef intable_type* intablePtr_type;
  typedef LookUpTable<KEY1,intablePtr_type> table_type;

public:

  /// Constructor
  LookUpTable2D(size_t table1size = 0,
  	size_t table2size = 0);

  /// Default destructor
  ~LookUpTable2D();

  /// Reserve memory
  /// @param size of the map to be set before starting inserting pairs in the  map
  /// @post  the memory corresponding to the given size will be reserved for
  ///        the future insertions.
  void reserve(size_t table1size,
         size_t table2size)
  {
    _table.reserve(table1size);
    _table2size = table2size;
  }

  /// Gets the interpolated VALUE of the supplied KEY's.
  /// @param key1 entry KEY1 to be interpolated
  /// @param key2 entry KEY2 to be interpolated
  VALUE get(const KEY1& key1, const KEY2& key2);

  /// Inserts a VALUE with a KEY as entry
  /// @param key1 the KEY1 to the value
  /// @param key2 the KEY2 to the value
  /// @param value the value to be entered
  void insert(const KEY1& key,
      const KEY2& key,
      const VALUE& value);

  /// Sorts the keys in the map
  void sortKeys();

  /// Outputs the table to the supplied stream
  /// @param outStream the stream to output the table
  void saveToStream(std::ostream& outStream);

  /// Inputs the table from the supplied stream
  /// @param inStream the stream to input the table from.
  void readFromStream(std::istream& inStream);

private:

  /// size of the inner table
  size_t _table2size;

  /// The actual map storage
  table_type _table;

}; // end of class LookUpTable2D

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "LookUpTable2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_LookUpTable2D_ci
