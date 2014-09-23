// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_ArrayT_hh
#define COOLFluiD_CudaEnv_ArrayT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CUDA/MacrosET.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
  
//////////////////////////////////////////////////////////////////////////////

#define ETVEC(a,b,c,d)   ArrayT<a<b<c,d>, ETuple<c,d,0> > >
#define ETMAT(a,b,c,d,e) ArrayT<a<b<c,d,e>, ETuple<c,d,e> > >
#define DATAPTR(a)       this->getData()->ptr()[a]

//////////////////////////////////////////////////////////////////////////////

/// Definition of a class ArrayT that implements an expression template technique
/// with statically defined size
/// @author Andrea Lani
template <typename DERIVED>
class ArrayT : public DERIVED {
public:
  typedef typename DERIVED::TUPLE::TYPE RET; 
  
  /// Default constructor
  HOST_DEVICE ArrayT(typename DERIVED::PTR* der) : DERIVED(der) {}
  
  /// Default destructor
  HOST_DEVICE virtual ~ArrayT() {}

  /// Overloading of the "[]" operator for assignment (writing).
  /// @param iElem index  
  HOST_DEVICE RET& operator[] (size_t iElem) {return DATAPTR(iElem);}

  /// Overloading of the "[]" operator for assignment (reading only).
  /// @param iElem index
  HOST_DEVICE RET operator[] (size_t iElem) const {return DATAPTR(iElem);}
  
  /// Accessor to individual entry
  /// @param iElem index
  HOST_DEVICE RET at (size_t iElem) const {return DATAPTR(iElem);}
  
  /// return the array size 
  HOST_DEVICE size_t size() const {return this->getData()->size();}
  
  /// @return maximum entry in the array
  HOST_DEVICE RET max() const 
  {
    RET emax = DATAPTR(0);
    VECLOOP(i, 1, size(), emax = (emax > DATAPTR(i)) ? emax : DATAPTR(i)); 
    return emax;
  }
  
  /// @return minimum entry in the array
  HOST_DEVICE RET min() const 
  {
    RET emin = DATAPTR(0);
    VECLOOP(i, 1, size(), emin = (emin < DATAPTR(i)) ? emin : DATAPTR(i));
    return emin;
  }
  
  /// @return sum of entry from start to end
  HOST_DEVICE RET sum() const 
  {
    RET esum = RET();
    VECLOOP(i, 0, size(), esum += DATAPTR(i)); 
    return esum;
  }

  /// @return sum of entry from start to end
  HOST_DEVICE RET partialSum(size_t start, size_t end) const 
  {
    RET esum = RET();
    VECLOOP(i, start, end, esum += DATAPTR(i)); 
    return esum;
  }  

  /// @return norm2 of entry from start to end
  HOST_DEVICE RET norm2() const 
  {
    return sqrt((double)(sqrNorm()));
  }
  
  /// @return norm1 of entry from start to end
  HOST_DEVICE RET norm1() const 
  {
    RET norm = abs(DATAPTR(0));
    for (size_t i = 1; i < size(); ++i)
      norm += abs(DATAPTR(i));
    return norm;
  } 
   
  /// @return normInf of entry from start to end
  HOST_DEVICE RET normInf() const 
  {
    RET norm = abs(DATAPTR(0));
    for (size_t i = 1; i < size(); ++i) {
      norm = (norm > abs(DATAPTR(i))) ? norm : abs(DATAPTR(i));
    }
    return norm;
  }  

  /// @return square of the norm2 of entry from start to end
  HOST_DEVICE RET sqrNorm() const 
  {
    RET sqsum = RET();
    VECLOOP(i, 0, size(), sqsum += DATAPTR(i)*DATAPTR(i));
    return sqsum;
  }   
  
protected:
  
  /// initialize the array with 0
  HOST_DEVICE void initialize(RET value) {VECLOOP(i, 0, size(), DATAPTR(i) = value);}
  
};
  
//////////////////////////////////////////////////////////////////////////////

 } // namespace CudaEnv

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
