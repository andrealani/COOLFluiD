// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ARRAYALLOCATOR_HH
#define ARRAYALLOCATOR_HH

///
/// ArrayAllocator
///
///  TODO:
///    - Use templates to add automagic padding?
///    - Implement removal of points (hint: free list inside the nodes)
///    - Make object-safe. ?? Is this needed?
///       --> Seems like Meshreader relies on it!! (strange)


#include "Common/BigAllocator.hh"
#include "Common/Obj_Helper.hh"
#include "Common/CFLog.hh"

namespace COOLFluiD {
    namespace Common {


template <class T, typename ALLOC = BigAllocator>
class ArrayAllocator
{
public:
  
  typedef CFuint IndexType;
  typedef T TYPE;
  
  /// Constructor
  ///   It allocates InitialSize elements of ElementSize bytes
  ArrayAllocator (const T & init, IndexType InitialSize,
    unsigned int ElementSize = sizeof (T));

  /// destructor
  ~ArrayAllocator ();
  
  ///   It allocates InitialSize elements of ElementSize bytes
  void initialize (const T & init, IndexType InitialSize,
                   unsigned int ElementSize = sizeof (T));
  
  ///   Clean up some local data built during initialization
  void free();
  
  /// Resize the array
  ///   !!! This can change MOVE the elements in memory!!!!
  ///   (but you're not supposed to apply the &-operator to
  ///   an element of this array)
  void Resize (IndexType NewSize);
  
  /// Grow the array by a number of elements
  /// (at least 1)
  /// The new number of elements is selected in function of
  /// resizing cost and granularity
  void grow ();
  
  /// Access to the array
  /// !!! Don't take the adress of an element!!!
  T & operator() (IndexType index);
  
  /// Access to the array
  /// !!! Don't take the adress of an element!!!
  const T & operator() (IndexType index) const;
   
  /// Access the address of the first entry
  T* ptr();
  
  /// Try to guess the index of the element pointed to
  IndexType PointerToIndex (const void * Ptr) const;
  
  /// For valarray compatibility
  size_t size () const;
  
  /// factor for which dividing the element size
  size_t sizeFactor() const {return 1;}
  
  /// For valarray compatibility
  void resize (size_t NewSize);
  
private:

  /// Disallow copy
  ArrayAllocator (const ArrayAllocator & A);

  /// Disallow assignment
  ArrayAllocator & operator = (const ArrayAllocator & A);
  
  /// This means the index is allowed from 0 .. getCurrSize()()-1
  IndexType getCurrSize() const;
  
private:
  ALLOC MemAlloc;
  IndexType CurrentSize;

  void * cf_restrict DataPtr;
  unsigned int ElementSize;
  
  bool CheckIndex(IndexType I) const;
  unsigned int ByteIndex(const IndexType I) const;
  T * PointerToIndex (const IndexType I) const;
};

///  Implementation
///***************************************************************************/

/// For valarray compatibility
template <class T, typename A>
inline size_t ArrayAllocator<T,A>::size () const {return getCurrSize();}

/// For valarray compatibility
template <class T, typename A>
inline void ArrayAllocator<T,A>::resize (size_t NewSize) {Resize (NewSize);}


template <class T, typename A>
inline unsigned int ArrayAllocator<T,A>::ByteIndex(const IndexType I) const
{
    cf_assert(I<CurrentSize);
    cf_assert(DataPtr != CFNULL);
    return (I*ElementSize);
}


template <class T, typename A>
inline T * ArrayAllocator<T,A>::PointerToIndex (const IndexType I) const
{
    return (T *) (static_cast<char *>(DataPtr) + ByteIndex(I));
}

template <class T, typename A>
inline T & ArrayAllocator<T,A>::operator() (IndexType index)
{
  CheckIndex(index);
  return *(T *) (static_cast<char *>(DataPtr)+ByteIndex(index));
}
      
template <class T, typename A>
inline const T & ArrayAllocator<T,A>::operator() (IndexType index) const
{
  CheckIndex(index);
  return *(T*)((T *)((char *) (DataPtr)+ByteIndex(index)));
}

template <class T, typename A>
inline T* ArrayAllocator<T,A>::ptr ()
{
  CheckIndex(0);
  return (T *) (static_cast<char *>(DataPtr)+ByteIndex(0));
}
      
template <class T, typename A>
typename ArrayAllocator<T,A>::IndexType ArrayAllocator<T,A>::PointerToIndex
(const void * Ptr) const
{
  // We have to cast to avoid pointer arithmetic corrections
  const char * BasePtr = reinterpret_cast<const char *>(DataPtr);
  const char * ElePtr = reinterpret_cast<const char *>(Ptr);
  
  if (ElePtr < BasePtr)
    return IndexType(-1);
  
  unsigned long Diff = ElePtr - BasePtr;
  if (Diff % ElementSize)
    return static_cast<typename ArrayAllocator<T,A>::IndexType>(-1);
  
  Diff /= ElementSize;
  if (Diff >= getCurrSize())
    return static_cast<IndexType>(-1);
  
  return IndexType(Diff);
}
      
template<class T, typename A>
inline bool ArrayAllocator<T,A>::CheckIndex(IndexType I) const
{
  cf_assert(I<getCurrSize());
  return true;
}      
      
template<class T, typename A>
void ArrayAllocator<T,A>::grow ()
{
  IndexType NewSize;
  const unsigned int Mingrow = 1;  // Expand at least this much elements
  unsigned int grow = MemAlloc.GetGranularity () / ElementSize;
  grow = (grow > Mingrow ? grow : Mingrow);
  
  cf_assert (grow != 0);
  
  if (MemAlloc.IsZeroCopy ()) {
    NewSize=CurrentSize+grow;
  }
  else {
    // If our current size is zero,
    // we grow at least one element and at most
    // up to the granularity
    NewSize=(CurrentSize ? 2*CurrentSize : grow);
  }
  cf_assert (NewSize > CurrentSize);
  Resize(NewSize);
}
      
/// Return the capacity (NUMBER OF ELEMENTS!)
template<class T, typename A>
inline
typename ArrayAllocator<T,A>::IndexType
ArrayAllocator<T,A>::getCurrSize () const
{
  return CurrentSize;
}

template<class T, typename A>
void ArrayAllocator<T,A>::Resize (IndexType count)
{
  IndexType OldSize = CurrentSize;
  
  // count is newsize, CurrentSize is current size
  if ((count < CurrentSize) &&
      (!Obj_Helper<T>::IsFundamental ()))
    {
      for (IndexType I = count; I<CurrentSize; I++)
	Obj_Helper<T>::destruct(PointerToIndex(I));
    }
  
  CFLogDebugMin( "Resizing ArrayAllocator " << (void*) this
      << " from " << OldSize << " to " << count << " ("
		 << count*ElementSize << " bytes)\n");
  
  //
  // !!! Resizing can move the POINTER !!!
    //
  MemAlloc.Resize (count*ElementSize);
  DataPtr=static_cast<void *> (MemAlloc.GetPtr ());
  cf_assert (!count || DataPtr);
  
  CurrentSize=count;
  
  //
  // For object safety
  //
  for (IndexType I = OldSize; I<CurrentSize; I++)
    Obj_Helper<T>::construct( PointerToIndex (I) );
}

template <class T, typename ALLOC>
ArrayAllocator<T,ALLOC>::ArrayAllocator (const T & Init,
					 IndexType InitSize,
					 unsigned int ESize)
  : CurrentSize (0), DataPtr(0), ElementSize(ESize)
{
  initialize(Init, InitSize, ESize);
}
      
template <class T, typename ALLOC>
void ArrayAllocator<T,ALLOC>::initialize (const T & Init,
					  IndexType InitSize,
					  unsigned int ESize)
{
  // reset the value of ElementSize
  ElementSize = ESize;

  // cf_assert (ElementSize >= sizeof (T));
  if (ElementSize >= sizeof (T)) {
    Resize (InitSize);
    cf_assert (getCurrSize () >= InitSize);

    // Init
    for (IndexType i=0; i < getCurrSize (); i++)
      Obj_Helper<T>::construct ( PointerToIndex(i), Init );

    CFLogDebugMin( "ArrayAllocator: init: size=" << InitSize
    << ", ElementSize=" << ElementSize << ", type-size="
    << sizeof (T) << "\n");
  }
}

template <class T, typename ALLOC>
void ArrayAllocator<T,ALLOC>::free()
{
  // This becomes an empty (hopefully optimized away)
  // loop in the case of primitive types
  for (IndexType I = 0; I < getCurrSize (); I++)
    Obj_Helper<T>::destruct ( PointerToIndex(I) );
}

template<class T, typename A>
ArrayAllocator<T,A>::~ArrayAllocator ()
{
  free();
}

} // namespace Common

} // COOLFluiD
#endif
