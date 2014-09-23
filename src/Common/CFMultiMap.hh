// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFMultiMap_hh
#define COOLFluiD_Common_CFMultiMap_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/NoSuchValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a "cheap" (from the point of view of the memory
/// requirements) map with MULTIPLE keys that gives the same
/// performance in searching of a std::multimap (the number
/// of comparisons is <= 2*log(size)+1. Since the high memory storage is
/// an issue with the standard <map>, this class offer a valid efficient
/// alternative to it.
/// Consider that once that you allocate a map, you can't get rid of the
/// allocated memory till the end of your run! This doesn't apply to the current
/// CFMultiMap that is just an Adapter of std::vector<std::pair<KEY,VALUE> >.
/// THIS IS REALLY AN IMPORTANT ADVANTAGE OF THIS MAP: since it is basically
/// a std::vector, YOU CAN RELEASE THE MEMORY when you go out of scope or you explicitly
/// delete a pointer to this object.
/// The class is parameterized with the types of the key and the type of the value
/// @pre this map can handle KEYs for which the operators "<" and ">" are available
/// @pre before using find() the CFMultiMap has to be sorted
/// @pre the situation in which this map must be used is a situation in which you first
///      insert ALL the pairs, then you sort (sortKeys() is a demanding operation
///      that must be applied only once!!) and then you start searching with find().
/// @author Andrea Lani
/// @author Tiago Quintino
template <class KEY, class VALUE>
class CFMultiMap {
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
  CFMultiMap(size_t maxSize = 0);

  /// Default destructor
  ~CFMultiMap();

  /// Reserve memory
  /// @param size of the map to be set before starting inserting pairs in the  map
  /// @post  the memory corresponding to the given size will be reserved for
  ///        the future insertions.
  void reserve(const size_t& maxSize)
  {
    _vectorMap.reserve(maxSize);
  }

  /// Clear the content of the map
  void clear()
  {
    std::vector<std::pair<KEY,VALUE> >().swap(_vectorMap);
  }

  /// Insert pair
  /// @param key   new key to be inserted
  /// @param value new value to be inserted, corresponding to the given key
  /// @post a new std::pair<KEY, VALUE> will be created and pushed back in the map
  /// @post the capacity of the map will increase only if the memory hasn't been
  ///       reserved properly at start up.
  void insert(const KEY& key, const VALUE& value);

  /// Find all the std::pair<KEY,VALUE> with the given key
  /// @param key      key to be looked-up
  /// @param isFound  flag telling if the value has been found
  /// @pre before using find() the CFMultiMap has to be sorted with sortKeys()
  /// @return a pair of iterators: the first iterator points to the first
  ///         pair having the requested key, the second iterator points to the
  ///         pair AFTER the last one sharing the same key. In fact, after the sort()
  ///         operation all the pairs having the same key are put one after the
  ///         other
 typedef typename std::vector<std::pair<KEY, VALUE> >::iterator MapIterator;
  std::pair<MapIterator, MapIterator> find(const KEY& key, bool& isFound);

  /// Get the number of pairs already inserted
  size_t getSize() const {return _vectorMap.size();}

  /// Get the number of pairs already inserted
  size_t size() const {return _vectorMap.size();}

  /// Sort all the pairs in the map by key
  /// @pre before using find() the CFMultiMap has to be sorted
  void sortKeys();

  /// Get the key in the current sequence corresponding to position i
  KEY getKey(CFuint i)
  {
    return _vectorMap[i].first;
  }

  /// Overloading of the operator"[]" for assignment
  /// @return _vectorMap[i].second
  VALUE& operator[] (const CFuint i)
  {
    cf_assert(i < size());
    return _vectorMap[i].second;
  }

  /// @return the iterator pointing at the first element
  MapIterator begin()
  {
    return _vectorMap.begin();
  }

  /// @return the iterator pointing at the last element
  MapIterator end()
  {
    return _vectorMap.end();
  }

private: // helper functions

  /// Error reporting in case of the key not being present
  void valueNotFound(const KEY& aKey) throw(Common::NoSuchValueException);

private: //nested classes

//////////////////////////////////////////////////////////////////////////////

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
  };

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
    /// @param p2   pair whose key is used as second term of comparison
    /// @return true  if(p1.first > p2.first)
    /// @post this is the second test to see if p1.first is == key during the 
    /// search with find()
    bool operator() (const std::pair<KEY,VALUE>& p1, 
		     const std::pair<KEY,VALUE>& p2) const
    {
      return (p1.first < p2.first);
    }
    
  };

private: //data

  /// Keeps track of the validity of the map
  /// if a key has been inserted and the map
  /// has not yet been sorted, then it this will be false
  bool _sorted;

  /// storage of the inserted data
  std::vector<std::pair<KEY, VALUE> > _vectorMap;

}; // end of class CFMultiMap

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFMultiMap.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFMultiMap_ci
