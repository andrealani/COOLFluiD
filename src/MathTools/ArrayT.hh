// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_ArrayT_hh
#define COOLFluiD_MathTools_ArrayT_hh

//////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <cassert>
#include "MathTools/MacrosET.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
  
//////////////////////////////////////////////////////////////////////////////

#define EETVEC(a,b,c,d)   ArrayT<a<b<c,d>, ETuple<c,d,0> > >
#define EETMAT(a,b,c,d,e) ArrayT<a<b<c,d,e>, ETuple<c,d,e> > >
#define EDATAPTR(a)       this->getData()->ptr()[a]

//////////////////////////////////////////////////////////////////////////////

/// Definition of a class ArrayT that implements an expression template technique
/// with statically defined size
/// @author Andrea Lani
template <typename DERIVED>
class ArrayT : public DERIVED {
public:
  typedef typename DERIVED::TUPLE::TYPE RET; 
  
  /// Default constructor
  HHOST_DEV ArrayT(typename DERIVED::PTR* der) : DERIVED(der) {}
  
  /// Default destructor
  HHOST_DEV virtual ~ArrayT() {}

  /// Overloading of the "[]" operator for assignment (writing).
  /// @param iElem index  
  HHOST_DEV RET& operator[] (size_t iElem) {assert(iElem < size()); return EDATAPTR(iElem);}

  /// Overloading of the "[]" operator for assignment (reading only).
  /// @param iElem index
  HHOST_DEV RET operator[] (size_t iElem) const {assert(iElem < size()); return EDATAPTR(iElem);}
  
  /// Accessor to individual entry
  /// @param iElem index
  HHOST_DEV RET at (size_t iElem) const {return EDATAPTR(iElem);}
  
  /// return the array size 
  HHOST_DEV size_t size() const {return this->getData()->size();}
  
  /// @return maximum entry in the array
  HHOST_DEV RET max() const 
  {
    RET emax = EDATAPTR(0);
    EVECLOOP(i, 1, size(), emax = std::max(emax, EDATAPTR(i))); 
    return emax;
  }
  
  /// @return minimum entry in the array
  HHOST_DEV RET min() const 
  {
    RET emin = EDATAPTR(0);
    EVECLOOP(i, 1, size(), emin = std::min(emin, EDATAPTR(i)));
    return emin;
  }
  
  /// @return sum of entry from start to end
  HHOST_DEV RET sum() const 
  {
    RET esum = RET();
    EVECLOOP(i, 0, size(), esum += EDATAPTR(i)); 
    return esum;
  }  
  
  /// @return sum of entry from start to end
  HHOST_DEV RET partialSum(size_t start, size_t end) const 
  {
    RET esum = RET();
    EVECLOOP(i, start, end, esum += EDATAPTR(i)); 
    return esum;
  }  

  /// @return norm2 of entry from start to end
  HHOST_DEV RET norm2() const 
  {
    return std::sqrt(sqrNorm());
  }
  
  /// @return norm1 of entry from start to end
  HHOST_DEV RET norm1() const 
  {
    assert(size() > 0);
    RET norm = std::abs(EDATAPTR(0));
    for (size_t i = 1; i < size(); ++i)
      norm += std::abs(EDATAPTR(i));
    return norm;
  } 
   
  /// @return normInf of entry from start to end
  HHOST_DEV RET normInf() const 
  {
    assert(size() > 0);
    RET norm = std::abs(EDATAPTR(0));
    for (size_t i = 1; i < size(); ++i)
      norm = std::max(norm, std::abs(EDATAPTR(i)));
    return norm;
  }  

  /// @return square of the norm2 of entry from start to end
  HHOST_DEV RET sqrNorm() const 
  {
    RET sqsum = RET();
    EVECLOOP(i, 0, size(), sqsum += EDATAPTR(i)*EDATAPTR(i));
    return sqsum;
  }   
  
protected:
  
  /// initialize the array with 0
  HHOST_DEV void initialize(RET value) {EVECLOOP(i, 0, size(), EDATAPTR(i) = value);}
  
};
  
//////////////////////////////////////////////////////////////////////////////

 } // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
