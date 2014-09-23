// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef GROWARRAY_HH
#define GROWARRAY_HH

#include "Common/ArrayAllocator.hh"
#include "Common/BigAllocator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common  {

//////////////////////////////////////////////////////////////////////////////


/// A primitive class that provides an dynamic resizeable array
/// that keeps track of the number of inserted elements.
/// This class is safe for non-primitive types (like classes)
/// but there is no requirement that the number of constructed objects is
/// the same as the number of objects in the GrowArray.
/// (It will be the capacity of the GrowArray)
    template <typename T, typename ALLOC = Common::BigAllocator>
    class GrowArray : public Common::ArrayAllocator<T, Common::BigAllocator>
    {
        typedef Common::ArrayAllocator<T, Common::BigAllocator> BaseClass;

        typedef typename BaseClass::IndexType IndexType;

      public:

        /// Construct a GrowArray of at least InitialSize capacity
        /// The number of elements will be InitialSize, and they will be
        /// initialised to T
        GrowArray ( const T & init, IndexType InitialSize,
                    size_t EleSize = sizeof ( T ) ) : BaseClass ( init, InitialSize, EleSize ), EleCount ( InitialSize )
        {
        }

        /// Return the number of objects in the array
        /// (This is can be changed by Resize() or increase())
        IndexType GetSize () const { return EleCount; }


        /// Return the number of objects that can be stored in the
        /// array without need for reallocation/growing
        IndexType GetCapacity () const { return BaseClass::GetSize (); }

        /// Resize the GrowArray
        /// This changes the number of objects in the array
        /// and also the capacity.
        /// The capacity will be at least NewSize, the number of
        /// objects will be exact NewSize
        void Resize ( IndexType NewSize )
        {
          BaseClass::Resize ( NewSize );
          EleCount = NewSize;
        }

        T & operator[] ( IndexType Index )
        {
          cf_assert ( Index < GetSize () );
          return BaseClass::operator() ( Index );
        }

        T & operator[] ( IndexType Index ) const
        {
          cf_assert ( Index < GetSize () );
          return BaseClass::operator() ( Index );
        }


        /// for compatibility
        /// TODO This does a default construct + one copy
        void push_back ( const T & Obj )
        {
          increase();
          ( *this ) [GetSize()-1] = Obj;
        }

        /// for compatibility with valarray
        size_t size () const { return GetSize (); }


        /// Make sure there is at least
        /// elecount capacity
        void reserve ( size_t Size )
        {
          if ( Size < GetSize () ) return ;
          BaseClass::Resize ( Size );
        }

        /// For compatibility with valarray
        void resize ( size_t NewSize ) { Resize ( NewSize ); }

        /// Remove all elements
        void clear () { EleCount = 0; }


        /// Increase the size by one and return the index of the new element
        /// If neccesairy, increase the capacity.
        IndexType increase ()
        {
          if ( GetSize () == GetCapacity () )
          {
            BaseClass::Grow ();
            cf_assert ( GetSize () < GetCapacity () );
          }
          EleCount++;
          return ( EleCount - 1 );
        }

      private:
        unsigned int EleCount;
    };


//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

#endif
