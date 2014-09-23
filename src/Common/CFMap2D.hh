// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFMap2D_hh
#define COOLFluiD_Common_CFMap2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NoSuchValueException.hh"
#include "Common/Trio.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class
/// @author Tiago Quintino
template <typename KEY1,
          typename KEY2,
          typename VALUE>
class CFMap2D {
public: // typedefs

  typedef typename std::vector< Trio<KEY1,KEY2,VALUE> >::iterator Iterator;

public: // methods

  /// Constructor
  CFMap2D(size_t maxSize = 0);

  /// Default destructor
  ~CFMap2D();

  /// Reserve memory
  /// @param size of the map to be set before starting inserting pairs in the  map
  /// @post  the memory corresponding to the given size will be reserved for
  ///        the future insertions.
  void reserve(const size_t& maxSize)
  {
    m_vectorMap.reserve(maxSize);
  }

  /// Insert pair
  /// @param key   new key to be inserted
  /// @param value new value to be inserted, corresponding to the given key
  /// @post a new std::pair<KEY, VALUE> will be created and pushed back in the map
  /// @post the capacity of the map will increase only if the memory hasn't been
  ///       reserved properly at start up.
  void insert(const KEY1& aKey1, const KEY2& aKey2, const VALUE& aValue);

  /// Find VALUE with the given KEY
  /// @param key  key to be looked-up
  /// @pre before using find() the CFMap2D has to be sorted with sortKeys()
  /// @post if the key is not in the map, a NoSuchValueException is thrown
  /// @exception NoSuchValueException is thrown if key not in map
  VALUE find(const KEY1& aKey1, const KEY2& aKey2);

  /// Find VALUE with the given KEY
  /// @param aKey  key to be looked-up
  /// @return the searched values corresponding to the inputed key
  /// @pre before using find() the CFMap3D has to be sorted with sortKeys()
  /// @post if the key is not in the map, a NoSuchValueException is thrown
  /// @exception NoSuchValueException is thrown if key not in map
  VALUE find(const std::pair<KEY1,KEY2>& aKey);

  /// Find the upper and lower pairs of KEY'S and VALUE'S bounding the
  /// supplied KEY.
  /// @param key  key to be looked-up
  /// @pre before using findBounds() the CFMap2D has to be sorted with sortKeys()
  /// @post if the key is not contained in the map what happens??
  /// @todo settle this question...
  std::pair<typename std::vector<Trio<KEY1,KEY2,VALUE> >::iterator,
             typename std::vector<Trio<KEY1,KEY2,VALUE> >::iterator>
  findBounds(const KEY1& aKey1, const KEY2& aKey2);

  /// Get the number of pairs already inserted
  size_t size() const
  {
    return m_vectorMap.size();
  }

  /// Overloading of the operator"[]" for assignment
  /// @return m_vectorMap[i].third
  VALUE& operator[] (const CFuint i)
  {
    cf_assert(i < size());
    return m_vectorMap[i].third;
  }

  /// Sort all the pairs in the map by key
  /// @pre before using find() the CFMap2D has to be sorted
  void sortKeys();

  /// Checks if a certain element is present
  /// @pre before using find() the CFMap3D has to be sorted
  bool isPresent(const Trio<KEY1,KEY2,VALUE>& value);

  /// @return the iterator pointing at the first element
  typename std::vector<Trio<KEY1, KEY2, VALUE> >::iterator begin()
  {
    return m_vectorMap.begin();
  }

  /// @return the iterator pointing at the last element
  typename std::vector<Trio<KEY1, KEY2, VALUE> >::iterator end()
  {
    return m_vectorMap.end();
  }

private: // helper functions

  /// Error reporting in case of the key not being present
  void valueNotFound(const KEY1& aKey1, const KEY2& aKey2);

private: // nested classes

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a functor object that is passed as an argument
  /// in the std::sort algo to compare two given pairs.
  /// @post the pairs are sorted in such a way that their KEYs are
  ///       listed in increasing order.
  /// @author Tiago Quintino
  class LessThan {
  public:

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param p1  pair to be used as first term of comparison
    /// @param p2  pair to be used as second term of comparison
    /// @return true  if(p1.first < p2.first)
    /// @return false if(p1.first >= p2.first)
    /// @post the pairs will be ordered according to the increasing order of
    ///       their keys
    /// @post sortKeys() uses these function to order the inserted pairs
    bool operator() (const Trio<KEY1,KEY2,VALUE>& p1,
                       const Trio<KEY1,KEY2,VALUE>& p2) const
    {
      bool ret;
      if (p1.first == p2.first) {
  ret = (p1.second < p2.second);
      }
      else {
  ret = (p1.first < p2.first);
      }
      return ret;
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param p1   pair whose key is used as first term of comparison
    /// @param keys  given key used as second term of comparison
    /// @post this is the first test to see if p1.first is == key during the
    ///       search with find()
    bool operator() (const Trio<KEY1,KEY2,VALUE>& p1, const std::pair<KEY1,KEY2>& keys) const
    {
      bool ret;
      if (p1.first == keys.first) {
  ret = (p1.second < keys.second);
      }
      else {
  ret = (p1.first < keys.first);
      }
      return ret;
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param keys  given key used as first term of comparison
    /// @param p1   pair whose key is used as second term of comparison
    /// @post this is the second test to see if p1.first is == key during the
    ///       search with find()
    bool operator() (const std::pair<KEY1,KEY2>& keys, const Trio<KEY1,KEY2,VALUE>& p1) const
    {
      bool ret;
      if (p1.first == keys.first) {
  ret = (p1.second > keys.second);
      }
      else {
  ret = (p1.first > keys.first);
      }
      return ret;
    }

  }; // end class LessThan

  /// This class represents a functor object that is passed as an argument
  /// in the std::sort algo to compare two given pairs.
  /// @post the pairs are sorted in such a way that their KEYs are
  ///       listed in increasing order.
  /// @author Tiago Quintino
  class Equal {
  public:

    /// Overloading of the operator() that makes this class acting as a functor
    /// @todo Missing documentation
    bool operator() (const Trio<KEY1,KEY2,VALUE>& p1,
                       const Trio<KEY1,KEY2,VALUE>& p2) const
    {
      return (p1.first == p2.first) && (p1.second == p2.second);
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @todo Missing documentation
    bool operator() (const Trio<KEY1,KEY2,VALUE>& p1,
           const std::pair<KEY1,KEY2>& keys) const
    {
      return (p1.first == keys.first) && (p1.second == keys.second);
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @todo Missing documentation
    bool operator() (const std::pair<KEY1,KEY2>& keys,
           const Trio<KEY1,KEY2,VALUE>& p1) const
    {
      return (p1.first == keys.first) && (p1.second == keys.second);
    }

  }; // end class Equal

//////////////////////////////////////////////////////////////////////////////

private: //data

  /// Keeps track of the validity of the map
  /// if a key has been inserted and the map
  /// has not yet been sorted, then it this will be false
  bool m_sorted;

  /// storage of the inserted data
  std::vector<Trio<KEY1,KEY2,VALUE> >  m_vectorMap;

}; // end of class CFMap2D

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFMap2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFMap2D_hh
