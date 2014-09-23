// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_CFMatSlice_hh
#define COOLFluiD_CudaEnv_CFMatSlice_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CUDA/ArrayT.hh"
#include "Common/CUDA/MatExprT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
  
  template <typename T, int N, int M> class CFMatSlice;

//////////////////////////////////////////////////////////////////////////////
  
/// Definition of a class CFMatSlice that implements an expression template technique
/// with statically defined size
/// @author Andrea Lani
template <typename T, int N = 0, int M = 0>
class CFMatSlice : public ETMAT(MatExprT,CFMatSlice,T,N,M) {
public:
  
  /// Default constructor
  HOST_DEVICE CFMatSlice(T* data, size_t totNbRows, size_t totNbCols) : 
    ETMAT(MatExprT,CFMatSlice,T,N,M)(this), 
    m_totNbRows(totNbRows), m_totNbCols(totNbCols), m_data(data) {}
    
  /// Copy constructor from other slice
  HOST_DEVICE CFMatSlice(const CFMatSlice<T,N,M>& in) : 
    ETMAT(MatExprT,CFMatSlice,T,N,M)(this), 
    m_totNbRows(in.m_totNbRows), m_totNbCols(in.m_totNbCols), m_data(in.m_data) {}
    
  /// Default destructor
  HOST_DEVICE ~CFMatSlice() {}
  
  /// Overloading of assignment operator(s)
#define MATSLICE_ASSIGN_OP(__op__) \
  template <typename EXPR>						\
    HOST_DEVICE const CFMatSlice<T,N,M>& operator __op__ (METYPEV(EXPR) expr) \
    {MATLOOP(i,0,N,j,0,M, (*this)(i COMA j) __op__ expr.at(i COMA j)); return *this;}
MATSLICE_ASSIGN_OP(=)
MATSLICE_ASSIGN_OP(+=)
MATSLICE_ASSIGN_OP(-=)
MATSLICE_ASSIGN_OP(*=)
MATSLICE_ASSIGN_OP(/=)
#undef MATSLICE_ASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define MATSLICE_ASSIGN_OP_CONST(__op__) \
  HOST_DEVICE const CFMatSlice<T,N,M>& operator __op__ (T expr) {VECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
MATSLICE_ASSIGN_OP_CONST(=)
MATSLICE_ASSIGN_OP_CONST(+=)
MATSLICE_ASSIGN_OP_CONST(-=)
MATSLICE_ASSIGN_OP_CONST(*=)
MATSLICE_ASSIGN_OP_CONST(/=)
#undef MATSLICE_ASSIGN_OP_CONST
  
  /// Overloading of the assignment operator "=" with CFVecSlice
  HOST_DEVICE const CFMatSlice<T,N,M>& operator= (const CFMatSlice<T,N,M>& other)
  {
    const size_t ns = size(); 
    for (size_t i = 0; i < ns; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// Overloading of the operator"()" for assignment.
  HOST_DEVICE T& operator() (size_t i, size_t j) {return m_data[i*m_totNbCols + j];}
 
  /// Overloading of the operator"()" (doesn't allow assignment)
  HOST_DEVICE T operator() (size_t i, size_t j) const {return m_data[i*m_totNbCols + j];}
  
  /// Overloading of the "[]" operator for assignment (writing).
  /// @param iElem index  
  HOST_DEVICE T& operator[] (size_t iElem) {return operator()(iElem/M, iElem%M);}
 
  /// Overloading of the "[]" operator for assignment (reading only).
  /// @param iElem index
  HOST_DEVICE T operator[] (size_t iElem) const {return operator()(iElem/M, iElem%M);}
  
  /// Accessor to individual entry
  /// @param iElem index
  HOST_DEVICE T at (size_t iElem) const {return operator()(iElem/M, iElem%M);}
  
  /// Accessor to individual entry
  /// @param i row index
  /// @param j column index    
  HOST_DEVICE T at (size_t i, size_t j) const {return operator()(i,j);}
  
  /// return the array size 
  HOST_DEVICE size_t size() const {return N*M;}
 
  /// return the number of rows 
  HOST_DEVICE size_t nbRows() const {return N;}
  
  /// return the number of columns 
  HOST_DEVICE size_t nbCols() const {return M;}
  
  /// copy content of another array
  template <typename ARRAY>
  HOST_DEVICE void copyFrom(const ARRAY& in){VECLOOP(i, 0, size(), m_data[i] = in[i]);}
  
  /// copy content of another array
  template <typename ARRAY>
  HOST_DEVICE void copyTo(ARRAY& in) {VECLOOP(i, 0, size(), in[i] = m_data[i]);}
    
  /// @return the raw data
  HOST_DEVICE T* ptr() {return m_data;}
  
private:
  
  /// total number of rows in the original matrix
  size_t m_totNbRows;
  
  /// total number of columns in the original matrix
  size_t m_totNbCols;
  
  /// array data
  T* m_data;
};
  
//////////////////////////////////////////////////////////////////////////////

/// Definition of a class CFMatSlice that implements an expression template technique
/// with dynamical size
/// @author Andrea Lani
template <typename T>
class CFMatSlice<T,0,0> : public ETMAT(MatExprT,CFMatSlice,T,0,0) {
public:
  
  // Constructor from preallocated memory
  HOST_DEVICE CFMatSlice(T* data, size_t totNbRows, size_t totNbCols, size_t ns, size_t ms) : 
    ETMAT(MatExprT,CFMatSlice,T,0,0)(this),  
    m_totNbRows(totNbRows), m_totNbCols(totNbCols), m_nrows(ns), m_ncols(ms), m_data(data) {}
    
  /// Copy constructor from other slice
  HOST_DEVICE CFMatSlice(const CFMatSlice<T,0,0>& in) : 
    ETMAT(MatExprT,CFMatSlice,T,0,0)(this),  
    m_totNbRows(in.m_totNbRows), m_totNbCols(in.m_totNbCols), m_nrows(in.m_nrows), m_ncols(in.m_ncols), m_data(in.m_data) {}
    
  /// Default destructor
  HOST_DEVICE ~CFMatSlice() {}
  
   /// Overloading of assignment operator(s)
#define MATSLICE0_ASSIGN_OP(__op__) \
  template <typename EXPR>						\
    HOST_DEVICE const CFMatSlice<T,0,0>& operator __op__ (METYPEV(EXPR) expr) \
    {MATLOOP(i,0,GETSIZER(NMAX(0,EXPR::SIZE1)),j,0,GETSIZEC(NMAX(0,EXPR::SIZE2)), (*this)(i COMA j) __op__ expr.at(i COMA j)); return *this;}
MATSLICE0_ASSIGN_OP(=)
MATSLICE0_ASSIGN_OP(+=)
MATSLICE0_ASSIGN_OP(-=)
MATSLICE0_ASSIGN_OP(*=)
MATSLICE0_ASSIGN_OP(/=)
#undef MATSLICE0_ASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define MATSLICE0_ASSIGN_OP_CONST(__op__) \
  HOST_DEVICE const CFMatSlice<T,0,0>& operator __op__ (T expr) {VECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
MATSLICE0_ASSIGN_OP_CONST(=)
MATSLICE0_ASSIGN_OP_CONST(+=)
MATSLICE0_ASSIGN_OP_CONST(-=)
MATSLICE0_ASSIGN_OP_CONST(*=)
MATSLICE0_ASSIGN_OP_CONST(/=)
#undef MATSLICE0_ASSIGN_OP_CONST

  /// Overloading of the assignment operator "=" with CFVecSlice
  HOST_DEVICE const CFMatSlice<T,0,0>& operator= (const CFMatSlice<T,0,0>& other)
  {
    const size_t ns = size(); 
    for (size_t i = 0; i < ns; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// Overloading of the operator"()" for assignment.
  HOST_DEVICE T& operator() (size_t i, size_t j) {return m_data[i*m_totNbCols + j];}

  /// Overloading of the operator"()" (doesn't allow assignment)
  HOST_DEVICE T operator() (size_t i, size_t j) const {return m_data[i*m_totNbCols + j];}
  
  /// Overloading of the "[]" operator for assignment (writing).
  /// @param iElem index  
  HOST_DEVICE T& operator[] (size_t iElem) {return operator()(iElem/m_ncols, iElem%m_ncols);}
 
  /// Overloading of the "[]" operator for assignment (reading only).
  /// @param iElem index
  HOST_DEVICE T operator[] (size_t iElem) const {return operator()(iElem/m_ncols, iElem%m_ncols);}
  
  /// Accessor to individual entry
  /// @param iElem index
  HOST_DEVICE T at (size_t iElem) const {return operator()(iElem/m_ncols, iElem%m_ncols);}
  
  /// Accessor to individual entry
  /// @param i row index
  /// @param j column index    
  HOST_DEVICE T at (size_t i, size_t j) const {return operator()(i,j);}
  
  /// return the array size 
  HOST_DEVICE size_t size() const {return m_nrows*m_ncols;}
 
  /// return the number of rows 
  HOST_DEVICE size_t nbRows() const {return m_nrows;}
  
  /// return the number of columns 
  HOST_DEVICE size_t nbCols() const {return m_ncols;}
  
  /// @return the raw data
  HOST_DEVICE T* ptr() {return m_data;}
  
private:
  
  /// total number of rows in the original matrix
  size_t m_totNbRows;
  
  /// total number of columns in the original matrix
  size_t m_totNbCols;
    
  /// number of matrix rows
  size_t m_nrows;
  
  /// number of matrix columns
  size_t m_ncols;
  
  /// array data
  T* m_data;
};
  
//////////////////////////////////////////////////////////////////////////////

  } // namespace CudaEnv

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
