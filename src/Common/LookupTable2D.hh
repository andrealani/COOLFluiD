// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_LookupTable2D_hh
#define COOLFluiD_Common_LookupTable2D_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a table where values can be looked up
/// and if not found, are interpolated with the nearest ones.
/// @todo Add trait class for interpolation algorithm
/// @author Andrea Lani
template <class KEY1, class KEY2, class VALUE>
class LookupTable2D {

public:

  /// Constructor
  LookupTable2D();

  /// Default destructor
  ~LookupTable2D();

  /// Initialize the table by reserving memory and
  /// creating the keys list
  /// @param key1   array containing all the key1s
  /// @param key2   array containing all the key2s
  /// @param valueStride stride for the value (that is
  ///                    logically an array of values)
  /// @pre this function MUST be called before using
  ///      insert() for the first time !!!
  void initialize(const std::vector<KEY1>& key1,
  	  const std::vector<KEY2>& key2,
  	  size_t valueStride);

  /// Gets the interpolated VALUE of the supplied KEY's.
  /// The interpolation strategy is based on the linear shape
  /// functions for a quadrilateral element
  /// @param key1 entry KEY1 to be interpolated
  /// @param key2 entry KEY2 to be interpolated
  VALUE get(const KEY1& key1,
    const KEY2& key2,
    const size_t iValue);

  /// Inserts a VALUE with a KEY as entry
  /// @param key1 the KEY1 to the value
  /// @param key2 the KEY2 to the value
  /// @param value the value to be entered
  void insert(const KEY1& key1,
      const KEY2& key2,
      const size_t iValue,
      const VALUE& value);

private:

  /// Get the look up index
  size_t getLookUpIndex(const KEY1& key1,
  const KEY2& key2,
  const size_t iValue)
  {
    size_t idx1 = _key1ToIdx.find(key1);
    size_t idx2 = _key2ToIdx.find(key2);
    return _valueStride*(idx1 + idx2*_key1ToIdx.size())
      + iValue;
  }

private:

  /// flag checking if the table is initialized
  bool _isInitialized;

  /// stride of the value array
  size_t _valueStride;

  /// Index map for the first key
  Common::CFMap<KEY1, size_t> _key1ToIdx;

  /// Index map for the second key
  Common::CFMap<KEY1, size_t> _key2ToIdx;

  /// storage of the values
  std::vector<VALUE> _values;

  /// storage of the temporary values
  std::vector<VALUE> _v;

}; // end of class LookupTable2D

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "LookupTable2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_LookupTable2D_ci
