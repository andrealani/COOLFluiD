// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFMap_hh
#define COOLFluiD_Common_CFMap_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/NoSuchValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a "cheap" (from the point of view of the memory
/// requirements) map with SINGLE key  that gives the same
/// performance in searching of a std::map (the number
/// of comparisons is <= 2*log(size)+1. Since the high memory storage is
/// an issue with the standard <map>, this class offer a valid efficient
/// alternative to it.
/// Consider that once that you allocate a map, you can't get rid of the
/// allocated memory till the end of your run! This doesn't apply to the current
/// CFMap that is just an Adapter of std::vector<std::pair<KEY,VALUE> >.
/// THIS IS REALLY AN IMPORTANT ADVANTAGE OF THIS MAP: since it is basically
/// a std::vector, YOU CAN RELEASE THE MEMORY when you go out of scope or you explicitly
/// delete a pointer to this object.
/// The class is parameterized with the types of the key and the type of the value
/// @pre this map can handle KEYs for which the operators "<" and ">" are available
/// @pre before using find() the CFMultiMap has to be sorted
/// @pre the situation in which this map must be used is a situation in which you first
///      insert ALL the pairs, then you sort (sortKeys() is a demanding operation
///      that must be applied only once!!) and then you start searching with find().
/// Consider using 2 seperate vectors (one for keys and one for elements)
/// It will be easier on the memory allocator, and it will provide for faster
/// access (multiple keys will be loaded in one cacheline when searching)
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename KEY, typename VALUE>
class CFMap {

public:

  /// Constructor
  /// @param maxSize  should be the maximum (foreseen) size of the map,
  ///                 allowing to efficiently allocate all the needed memory,
  ///                 without having to do that on the fly while inserting new
  ///                 pairs in the map.
  /// @post           The default maximum size is set to 0. The function reserve()
  ///                 allow you to set the size not necessarly at construction time.
  ///                 You can start inserting without specifying a size neither
  ///                 in the constructor, nor with reserve(), but this will be
  ///                 much more inefficient.
  /// @post           if the maximum size is supplied, the maxSize
  ///                 of the map will be set to maxSize.
  CFMap(size_t maxSize = 0);

  /// Default destructor
  ~CFMap();

  /// Reserve memory
  /// @param size of the map to be set before starting inserting pairs in the  map
  /// @post  the memory corresponding to the given size will be reserved for
  ///        the future insertions.
  void reserve(size_t maxSize)
  {
    m_vectorMap.reserve(maxSize);
  }

  /// Insert pair
  /// @param key   new key to be inserted
  /// @param value new value to be inserted, corresponding to the given key
  /// @post a new std::pair<KEY, VALUE> will be created and pushed back in the map
  /// @post the capacity of the map will increase only if the memory hasn't been
  ///       reserved properly at start up.
  void insert(const KEY& aKey, const VALUE& aValue);

  /// Find VALUE with the given KEY
  /// @param key  key to be looked-up
  /// @pre before using find() the CFMap has to be sorted with sortKeys()
  /// @post if the key is not in the map, a NoSuchValueException is thrown
  /// @return copy of found VALUE
  /// @exception NoSuchValueException is thrown if key not in map
  VALUE find(const KEY& aKey);

  /// Find all the std::pair<KEY,VALUE> with the given key
  /// @param key      key to be looked-up
  /// @param isFound  flag telling if the value has been found
  /// @pre before using find() the CFMultiMap has to be sorted with sortKeys()
  /// @return copy of found VALUE
  VALUE find(const KEY& key, bool& isFound);

  /// Check if the given KEY is associated to a VALUE
  /// @param key  key to be looked-up
  /// @pre before using exists() the CFMap has to be sorted with sortKeys()
  /// @return flag to know if key exists
  bool exists(const KEY& aKey);

  /// Find the upper and lower pairs of KEY'S and VALUE'S bounding the
  /// supplied KEY.
  /// @param key  key to be looked-up
  /// @pre before using findBounds() the CFMap has to be sorted with sortKeys()
  /// @post if the key is not contained in the map what happens??
  /// @todo settle this question...
  std::pair<typename std::vector<typename std::pair<KEY,VALUE> >::iterator,
             typename std::vector<typename std::pair<KEY,VALUE> >::iterator>
  findBounds(const KEY& key);

  /// Print the content of the map
  void print();

  /// Clear the content of the map
  void clear()
  {
    std::vector<std::pair<KEY,VALUE> >().swap(m_vectorMap);
  } 
  
  /// Clear the content of the map while keeping its size
  void clearContent() 
  {
    if (size() > 0 ) {m_vectorMap.clear();}
    cf_assert(size() == 0);
  }
  
  /// Get the number of pairs already inserted
  size_t size() const
  {
    return m_vectorMap.size();
  }

  /// Overloading of the operator"[]" for assignment
  /// @return m_vectorMap[i].second
  VALUE& operator[] (const CFuint i)
  {
    cf_assert(i < size());
    return m_vectorMap[i].second;
  }

  /// Sort all the pairs in the map by key
  /// @pre before using find() the CFMap has to be sorted
  void sortKeys();

  /// @return the iterator pointing at the first element
  typename std::vector<std::pair<KEY, VALUE> >::iterator begin()
  {
    return m_vectorMap.begin();
  }

  /// @return the iterator pointing at the last element
  typename std::vector<std::pair<KEY, VALUE> >::iterator end()
  {
    return m_vectorMap.end();
  }

  /// replace the value of the key: Key, with the value: newValue
  /// @autor Alessandro Sanna
  void replaceValue(KEY Key, VALUE newValue)
  {
    if(!m_sorted) { sortKeys(); }
    VALUE oldValue = find(Key);
    std::pair<KEY, VALUE > oldPair(Key, oldValue);
    std::pair<KEY, VALUE > newPair(Key, newValue);
    replace(m_vectorMap.begin(), m_vectorMap.end(), oldPair, newPair);
  }

  /// get the key corresponding to position i
  /// @autor Alessandro Sanna
  KEY getKey(const CFuint i)
  {
    cf_assert(i < size());
    if(!m_sorted) { sortKeys(); }
    return m_vectorMap[i].first;
  }


private: // helper functions

  /// Error reporting in case of the key not being present
  void valueNotFound(const KEY& aKey);

private: // nested classes

  /// This class represents a functor object that is passed as an argument
  /// in the std::sort algo to compare two given pairs.
  /// @post the pairs are sorted in such a way that their KEYs are
  ///       listed in increasing order.
  /// @author Andrea Lani
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
    bool operator() (const std::pair<KEY,VALUE>& p1,
                       const std::pair<KEY,VALUE>& p2) const
    {
      return (p1.first < p2.first) ? true : false;
    }

  }; // end class LessThan

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a functor object that is used to compare a
  /// given key with the key in a std::pair. This functor is passed as
  /// an argument in the std::equal_range function in order to find all
  /// the pairs containing the specified key.
  /// @author Andrea Lani
  class Compare {
  public:

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param p1   pair whose key is used as first term of comparison
    /// @param key  given key used as second term of comparison
    /// @return true  if(p1.first < key)
    /// @return false if(p1.first >= key)
    /// @post this is the first test to see if p1.first is == key during the
    ///       search with find()
    bool operator() (const std::pair<KEY,VALUE>& p1, const KEY& key) const
    {
      return (p1.first < key) ? true : false;
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param key  given key used as first term of comparison
    /// @param p1   pair whose key is used as second term of comparison
    /// @return true  if(p1.first > key)
    /// @return false if(p1.first <= key)
    /// @post this is the second test to see if p1.first is == key during the
    ///       search with find()
    bool operator() (const KEY& key, const std::pair<KEY,VALUE>& p1) const
    {
      return (p1.first > key) ? true : false;
    }

  }; // class Compare

private: //data

  /// Keeps track of the validity of the map
  /// if a key has been inserted and the map
  /// has not yet been sorted, then it this will be false
  bool m_sorted;

  /// storage of the inserted data
  std::vector<std::pair<KEY,VALUE> >  m_vectorMap;

}; // end of class CFMap

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFMap.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFMap_hh
