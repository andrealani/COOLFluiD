// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_CFVecSlice_hh
#define COOLFluiD_CudaEnv_CFVecSlice_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CUDA/ExprT.hh"

// #ifndef CXX_NEEDS_FRIEND_TMPL_DECL
// #define CXX_NEEDS_FRIEND_TMPL_DECL
// #endif

#include "Common/Compatibility.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
  
  template  <typename T, int N> class CFVecSlice;
    
//////////////////////////////////////////////////////////////////////////////
  
/// Definition of a class CFVecSlice that implements an expression template technique
/// with statically defined size
/// @author Andrea Lani
template <typename T, int N = 0>
class CFVecSlice : public ETVEC(ExprT,CFVecSlice,T,N) {
public:
  
  /// Default constructor
  HOST_DEVICE CFVecSlice(T* data) : ETVEC(ExprT,CFVecSlice,T,N)(this), m_data(data) {}
    
  /// Copy constructor from other slice
  HOST_DEVICE CFVecSlice(const CFVecSlice<T,N>& in) : ETVEC(ExprT,CFVecSlice,T,N)(this), m_data(in.m_data) {}
      
  /// Default destructor
  HOST_DEVICE ~CFVecSlice() {}
  
  /// Overloading of assignment operator(s)
#define VECSLICE_ASSIGN_OP(__op__) \
  template <typename EXPR>						\
  HOST_DEVICE const CFVecSlice<T,N>& operator __op__ (ETYPEV(EXPR) expr) \
    {VECLOOP(i, 0, N, m_data[i] __op__ expr.at(i)); return *this;}
VECSLICE_ASSIGN_OP(=)
VECSLICE_ASSIGN_OP(+=)
VECSLICE_ASSIGN_OP(-=)
VECSLICE_ASSIGN_OP(*=)
VECSLICE_ASSIGN_OP(/=)
#undef VECSLICE_ASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define VECSLICE_ASSIGN_OP_CONST(__op__) \
  HOST_DEVICE const CFVecSlice<T,N>& operator __op__ (T expr) \
  {VECLOOP(i, 0, N, m_data[i] __op__ expr); return *this;}
VECSLICE_ASSIGN_OP_CONST(=)
VECSLICE_ASSIGN_OP_CONST(+=)
VECSLICE_ASSIGN_OP_CONST(-=)
VECSLICE_ASSIGN_OP_CONST(*=)
VECSLICE_ASSIGN_OP_CONST(/=)
#undef VECSLICE_ASSIGN_OP_CONST

  /// Overloading of the assignment operator "=" with CFVecSlice
  HOST_DEVICE const CFVecSlice<T,N>& operator= (const CFVecSlice<T,N>& other)
  {
    for (size_t i = 0; i < N; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// return the array size 
  HOST_DEVICE size_t size() const {return N;}

  /// copy content of another array
  template <typename ARRAY>
  HOST_DEVICE void copyFrom(const ARRAY& in){VECLOOP(i, 0, N, m_data[i] = in[i]);}
  
  /// copy content of another array
  template <typename ARRAY>
  HOST_DEVICE void copyTo(ARRAY& in) {VECLOOP(i, 0, N, in[i] = m_data[i]);}
    
  /// @return the raw data
  HOST_DEVICE T* ptr() {return m_data;}
  
private:
  
  /// array data
  T* m_data;
};
  
//////////////////////////////////////////////////////////////////////////////

/// Definition of a class CFVecSlice that implements an expression template technique
/// with dynamical size
/// @author Andrea Lani
template <typename T>
class CFVecSlice<T,0> : public ETVEC(ExprT,CFVecSlice,T,0) {
public:
  
  // Constructor from preallocated memory
  HOST_DEVICE CFVecSlice(T* data, size_t ns) : 
    ETVEC(ExprT,CFVecSlice,T,0)(this), m_size(ns), m_data(data) {}
    
  /// Copy constructor from other slice
  HOST_DEVICE CFVecSlice(const CFVecSlice<T,0>& in) : 
    ETVEC(ExprT,CFVecSlice,T,0)(this), m_size(in.size()), m_data(in.m_data) {}
    
  /// Default destructor
  HOST_DEVICE ~CFVecSlice() {}
  
   /// Overloading of assignment operator(s)
#define VECSLICE0_ASSIGN_OP(__op__) \
  template <typename EXPR>						\
  HOST_DEVICE const CFVecSlice<T,0>& operator __op__ (ETYPEV(EXPR) expr) \
    {VECLOOP(i, 0, GETSIZE1(NMAX(0,EXPR::SIZE1)), m_data[i] __op__ expr.at(i)); return *this;}
VECSLICE0_ASSIGN_OP(=)
VECSLICE0_ASSIGN_OP(+=)
VECSLICE0_ASSIGN_OP(-=)
VECSLICE0_ASSIGN_OP(*=)
VECSLICE0_ASSIGN_OP(/=)
#undef VECSLICE0_ASSIGN_OP

/// Overloading of assignment operator(s) with constants
#define VECSLICE0_ASSIGN_OP_CONST(__op__) \
  HOST_DEVICE const CFVecSlice<T,0>& operator __op__ (T expr) \
  {VECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
VECSLICE0_ASSIGN_OP_CONST(=)
VECSLICE0_ASSIGN_OP_CONST(+=)
VECSLICE0_ASSIGN_OP_CONST(-=)
VECSLICE0_ASSIGN_OP_CONST(*=)
VECSLICE0_ASSIGN_OP_CONST(/=)
#undef VECSLICE0_ASSIGN_OP_CONST
 
  /// Overloading of the assignment operator "=" with CFVecSlice
  HOST_DEVICE const CFVecSlice<T,0>& operator= (const CFVecSlice<T,0>& other)
  {
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// return the array size 
  HOST_DEVICE size_t size() const {return m_size;}
 
  /// @return the raw data
  HOST_DEVICE T* ptr() {return m_data;}
    
private:
  
  /// array size
  size_t m_size;
  
  /// array data
  T* m_data;
};
  
//////////////////////////////////////////////////////////////////////////////
    
  } // namespace CudaEnv

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
