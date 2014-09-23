// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFMap3D_hh
#define COOLFluiD_Common_CFMap3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NoSuchValueException.hh"
#include "Common/StringOps.hh"
#include "Common/Quartet.hh"
#include "Common/Trio.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class implements an associative container where the key is triple and
/// the corresponding value has multeplicity >= 1
/// @author Tiago Quintino
/// @author Andrea Lani
template <typename KEY1,
          typename KEY2,
          typename KEY3,
          typename VALUE>
class CFMap3D {
public: // typedefs

  typedef typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator MapIterator;

public: // methods

  /// Constructor
  CFMap3D(size_t maxSize = 0);

  /// Default destructor
  ~CFMap3D();

  /// Reserve memory
  /// @param size of the map to be set before starting inserting pairs in the  map
  /// @post  the memory corresponding to the given size will be reserved for
  ///        the future insertions.
  void reserve(const size_t& maxSize) {  m_vectorMap.reserve(maxSize); }

  /// Clear the content of the map
  void clear() {  std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >().swap(m_vectorMap); }

  /// Insert a value with corresponding keys
  /// @param aKey1 1st key
  /// @param aKey2 2nd key
  /// @param aKey3 3rd key
  /// @param aValue new value to be inserted, corresponding to the given keys
  /// @post the capacity of the map will increase only if the memory hasn't been
  ///       reserved properly at start up.
  void insert(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3, const VALUE& aValue);

  /// Insert a value with corresponding keys
  /// @param aKey a Trio representing the 3 keys
  /// @param aValue new value to be inserted, corresponding to the given keys
  /// @post the capacity of the map will increase only if the memory hasn't been
  ///       reserved properly at start up.
  void insert(const Trio<KEY1,KEY2,KEY3>& aKey, const VALUE& aValue);

  /// Find VALUE with the given KEY
  /// @param aKey1 1st key to locate
  /// @param aKey2 2nd key to locate
  /// @param aKey3 3rd key to locate
  /// @return the searched values corresponding to the inputed keys
  /// @pre before using find() the CFMap3D has to be sorted with sortKeys()
  /// @post if the key is not in the map, a NoSuchValueException is thrown
  /// @exception NoSuchValueException is thrown if key not in map
  VALUE find(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3);

  /// Find all the std::pair<KEY,VALUE> with the given key
  /// @param key      key to be looked-up
  /// @param isFound  flag telling if the value has been found
  /// @pre before using find() the CFMultiMap has to be sorted with sortKeys()
  /// @return a pair of iterators: the first iterator points to the first
  ///         pair having the requested key, the second iterator points to the
  ///         pair AFTER the last one sharing the same key. In fact, after the sort()
  ///         operation all the pairs having the same key are put one after the
  ///         other
  std::pair<MapIterator, MapIterator> find(const KEY1& aKey1,
  				   const KEY2& aKey2,
  				   const KEY3& aKey3,
  				   bool& isFound);

  /// Find VALUE with the given KEY
  /// @param aKey  key to be looked-up
  /// @return the searched values corresponding to the inputed key
  /// @pre before using find() the CFMap3D has to be sorted with sortKeys()
  /// @post if the key is not in the map, a NoSuchValueException is thrown
  /// @exception NoSuchValueException is thrown if key not in map
  VALUE find(const Trio<KEY1,KEY2,KEY3>& aKey);

  /// Find the upper and lower pairs of KEY'S and VALUE'S bounding the
  /// supplied KEY.
  /// @param key  key to be looked-up
  /// @pre before using findBounds() the CFMap3D has to be sorted with sortKeys()
  /// @post if the key is not contained in the map what happens??
  /// @todo settle this question...
  std::pair<typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator,
             typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator>
  findBounds(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3);

  /// Get the number of pairs already inserted
  size_t size() const { return m_vectorMap.size(); }

  /// Prints all the elements of the map to the ostream
  /// @param stream he stream to output to
  std::ostream& print(std::ostream& out)
  {
    for(CFuint i = 0; i < m_vectorMap.size(); ++i) {
      out << m_vectorMap[i] << "\n";
    }
    return out;
  }

  /// Overloading of the operator"[]" for assignment
  /// @return m_vectorMap[i].fourth
  VALUE& operator[] (const CFuint i)
  {
    cf_assert(i < size());
    return m_vectorMap[i].fourth;
  }

  /// Gets the i entry
  Quartet<KEY1,KEY2,KEY3,VALUE>& getEntry(const CFuint i)
  {
    cf_assert(i < size());
    return m_vectorMap[i];
  }

  /// Sort all the pairs in the map by key
  /// @pre before using find() the CFMap3D has to be sorted
  void sortKeys();

  /// Checks if an element with a certain key has been inserted
  bool exists (const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3);

  /// @return the iterator pointing at the first element
  typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator begin() { return m_vectorMap.begin(); }

  /// @return the iterator pointing at the last element
  typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator end() { return m_vectorMap.end(); }

private: // helper functions

  /// Error reporting in case of the key not being present
  void valueNotFound(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3);

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
    /// @return true  if p1.first less p2.first
    /// @return false if(p1.first >= p2.first)
    /// @post the pairs will be ordered according to the increasing order of
    ///       their keys
    /// @post sortKeys() uses these function to order the inserted pairs
    bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1,
                       const Quartet<KEY1,KEY2,KEY3,VALUE>& p2) const
    {
      bool ret;
      if (p1.first == p2.first) {
        if (p1.second == p2.second) {
          ret = (p1.third < p2.third);
        }
        else {
          ret = (p1.second < p2.second);
        }
      }
      else {
        ret = (p1.first < p2.first);
      }
      return ret;
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param p1   pair whose key is used as first term of comparison
    /// @param keys  given key used as second term of comparison
    /// @post this is the first test to see if p1.first is == key during the search with find()
    bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1, const Trio<KEY1,KEY2,KEY3>& keys) const
    {
      bool ret;
      if (p1.first == keys.first) {
        if (p1.second == keys.second) {
          ret = (p1.third < keys.third);
        }
        else {
          ret = (p1.second < keys.second);
        }
      }
      else {
        ret = (p1.first < keys.first);
      }
      return ret;
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @param keys  given key used as first term of comparison
    /// @param p1   pair whose key is used as second term of comparison
    /// @post this is the second test to see if p1.first is == key during the search with find()
    bool operator() (const Trio<KEY1,KEY2,KEY3>& keys, const Quartet<KEY1,KEY2,KEY3,VALUE>& p1) const
    {
      bool ret;
      if (p1.first == keys.first) {
        if (p1.second == keys.second) {
          ret = (p1.third > keys.third);
        }
        else {
          ret = (p1.second > keys.second);
        }
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
    bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1,
                     const Quartet<KEY1,KEY2,KEY3,VALUE>& p2) const
    {
      return (p1.first == p2.first) && (p1.second == p2.second) && (p1.third == p2.third);
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @todo Missing documentation
    bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1,
                     const Trio<KEY1,KEY2,KEY3>& keys) const
    {
      return operator()(keys,p1);
    }

    /// Overloading of the operator() that makes this class acting as a functor
    /// @todo Missing documentation
    bool operator() (const Trio<KEY1,KEY2,KEY3>& keys,
                     const Quartet<KEY1,KEY2,KEY3,VALUE>& p1) const
    {
      return (p1.first == keys.first) && (p1.second == keys.second) && (p1.third == keys.third);
    }

  }; // end class Equal

  struct EqualComp
  {
    const Trio<KEY1,KEY2,KEY3>& m_key;
    EqualComp (const Trio<KEY1,KEY2,KEY3>& aKey) : m_key (aKey) {}
    bool operator() ( const Quartet<KEY1,KEY2,KEY3,VALUE>& obj )
    {
      return (m_key.first == obj.first) && (m_key.second == obj.second) && (m_key.third == obj.third);
    };
  };


//////////////////////////////////////////////////////////////////////////////

private: //data

  /// Keeps track of the validity of the map
  /// if a key has been inserted and the map
  /// has not yet been sorted, then it this will be false
  bool m_sorted;

  /// storage of the inserted data
  std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >  m_vectorMap;

}; // end of class CFMap3D

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
 CFMap3D<KEY1, KEY2, KEY3, VALUE>::CFMap3D(size_t maxSize) :
   m_sorted(false)
{
  if(maxSize != 0) {
    m_vectorMap.reserve(maxSize);
  }
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
CFMap3D<KEY1, KEY2, KEY3, VALUE>::~CFMap3D()
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline void CFMap3D<KEY1, KEY2, KEY3, VALUE>::insert(const KEY1& aKey1,
                                        const KEY2& aKey2,
                                        const KEY3& aKey3,
                                        const VALUE& aValue)
{
  m_sorted = false;
  m_vectorMap.push_back(Quartet<KEY1, KEY2, KEY3, VALUE>(aKey1, aKey2, aKey3, aValue));
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline void CFMap3D<KEY1, KEY2, KEY3, VALUE>::insert(const Trio<KEY1,KEY2,KEY3>& aKey,
                                              const VALUE& aValue)
{
  insert(aKey.first, aKey.second, aKey.third, aValue);
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline std::pair<
  typename std::vector<Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator,
  typename std::vector<Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator>
CFMap3D<KEY1, KEY2, KEY3, VALUE>::findBounds(const KEY1& aKey1,
                                             const KEY2& aKey2,
                                             const KEY3& aKey3)
{
  if(m_vectorMap.empty()) {  valueNotFound(aKey1,aKey2,aKey3); }

  if(!m_sorted) { sortKeys(); }

  MapIterator itr =
   std::lower_bound(m_vectorMap.begin(),
                    m_vectorMap.end(),
                    make_Trio(aKey1,aKey2,aKey3),
                    LessThan());

  if((itr->first != aKey1) || (itr->second != aKey2) || (itr->third != aKey3)) {
    valueNotFound(aKey1,aKey2,aKey3);
  }

  Equal eq;
  Trio<KEY1,KEY2,KEY3> key(aKey1,aKey2,aKey3);
  MapIterator before = itr;
  for( ; itr != m_vectorMap.end(); ++itr) {
    if (!eq(*itr,key)) break;
  }

  return std::make_pair(before,itr);
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline void CFMap3D<KEY1, KEY2, KEY3, VALUE>::valueNotFound(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
{
   std::string msg = "CFMap3D: KEYS not found: ";
   msg += Common::StringOps::to_str(aKey1);
   msg += " ";
   msg += Common::StringOps::to_str(aKey2);
   msg += " ";
   msg += Common::StringOps::to_str(aKey3);
   throw Common::NoSuchValueException (FromHere(),msg);
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline VALUE CFMap3D<KEY1, KEY2, KEY3, VALUE>::find(const Trio<KEY1,KEY2,KEY3>& aKey)
{
  return find(aKey.first, aKey.second, aKey.third);
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline bool CFMap3D<KEY1, KEY2, KEY3, VALUE>::exists(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
{
  if(!m_sorted) { sortKeys(); }

  MapIterator itr =
   std::find_if ( m_vectorMap.begin(),
                  m_vectorMap.end(),
                  EqualComp ( make_Trio(aKey1,aKey2,aKey3 ) ) );

  return itr != m_vectorMap.end();
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline VALUE CFMap3D<KEY1, KEY2, KEY3, VALUE>::find(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
{
  if(m_vectorMap.empty()) {
    valueNotFound(aKey1,aKey2,aKey3);
  }

  if(!m_sorted) { sortKeys(); }

  MapIterator itr =
   std::lower_bound(m_vectorMap.begin(),
                    m_vectorMap.end(),
                    make_Trio(aKey1,aKey2,aKey3),
                    LessThan());

 if((itr->first != aKey1) || (itr->second != aKey2) || (itr->third != aKey3)) {
    valueNotFound(aKey1,aKey2,aKey3);
 }

 return itr->fourth;
}

//////////////////////////////////////////////////////////////////////////////

template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline std::pair<typename std::vector<Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator,
		 typename std::vector<Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator>
CFMap3D<KEY1, KEY2, KEY3, VALUE>::find(const KEY1& aKey1, const KEY2& aKey2,
				       const KEY3& aKey3, bool& isFound)
{
  if(!m_sorted) { sortKeys(); }
  
  std::pair<MapIterator,MapIterator> result;
  // AL: note that the map can have 0 size (in parallel programming it
  // can easily happen!)
  // this fix is needed because if you apply an equal_range to a vector
  // with null size you get crazy stuff ...
  // this MUST NOT be changed !!!
  if (m_vectorMap.empty()) {
    isFound = false;
    return result;
  }
  
  result=std::equal_range(m_vectorMap.begin(),
			  m_vectorMap.end(),
			  make_Trio(aKey1,aKey2,aKey3),
			  LessThan());
  
  isFound = !((result.first->first != aKey1) ||
	      (result.first->second != aKey2) ||
	      (result.first->third != aKey3));
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////
    
template<typename KEY1, typename KEY2, typename KEY3, typename VALUE>
inline void CFMap3D<KEY1, KEY2, KEY3, VALUE>::sortKeys()
{
  std::sort(m_vectorMap.begin(),m_vectorMap.end(), LessThan());
  m_sorted = true;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFMap3D_hh
