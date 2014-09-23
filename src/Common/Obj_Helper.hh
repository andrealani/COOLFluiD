// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_OBJ_HELPER_HH
#define COOLFluiD_Common_OBJ_HELPER_HH

#include <cstring>

#include "Common/CppTypeInfo.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common  {

//////////////////////////////////////////////////////////////////////////////

    /// Forward declaration
    template <class Type> struct Obj_Helper;

    /// Default: assume objects need to be constructed
    template <class Type, bool Simple>
    struct Obj_Helper_Helper
    {
    private:
      inline static void construct ( Type * memptr ) { new ( memptr ) Type(); }

      inline static void construct ( Type * memptr, const Type & copy ) { new ( memptr ) Type ( copy );     }

      inline static void destruct ( Type * memptr ) { memptr->~Type (); }

      friend struct Obj_Helper<Type>;
    };


    /// Optimized version for simple types
    template <class Type>
    struct  Obj_Helper_Helper<Type,true>
    {
    private:
      inline static void construct ( Type * memptr, const Type & copy ) { *memptr = copy; }

      inline static void construct ( Type * memptr )
      {
#     ifndef NDEBUG
        memset ( memptr, 0, sizeof ( Type ) );
#     endif
      }

      inline static void destruct ( Type * memptr ) { /* Do nothing */  }

      friend struct Obj_Helper<Type>;
    };

    /// Helper class for turning a raw memory area into a pool
    /// of objects (or in the case of simple types initialiase them fast
    /// Easy to use interface
    template <class Type>
    struct Obj_Helper
    {
      inline static void destruct ( Type * memptr )
      {
        Obj_Helper_Helper<Type, is_fundamental<Type>::type>
        ::destruct ( memptr );
      }

      inline static void construct ( Type * memptr )
      {
        Obj_Helper_Helper<Type, is_fundamental<Type>::type>
        ::construct ( memptr );
      }

      inline static void construct ( Type * memptr, const Type & init )
      {
        Obj_Helper_Helper<Type, is_fundamental<Type>::type>
        ::construct ( memptr, init );
      }

      inline static bool IsFundamental ()
      {
        return is_fundamental<Type>::type;
      };
    };

//////////////////////////////////////////////////////////////////////////////

  }  // Common
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_OBJ_HELPER_HH
