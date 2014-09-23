// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_CFMatSlice_hh
#define COOLFluiD_MathTools_CFMatSlice_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/ArrayT.hh"
#include "MathTools/MatExprT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
  
  template <typename T, int N, int M> class CFMatSlice;

  template <typename T, int N, int M> std::ostream& operator<< 
    (std::ostream& out, const CFMatSlice<T,N,M>& A);

  template <typename T, int N, int M> std::istream& operator>> 
    (std::istream& in, CFMatSlice<T,N,M>& A);
  
//////////////////////////////////////////////////////////////////////////////
  
/// Definition of a class CFMatSlice that implements an expression template technique
/// with statically defined size
/// @author Andrea Lani
template <typename T, int N = 0, int M = 0>
class CFMatSlice : public EETMAT(MatExprT,CFMatSlice,T,N,M) {
public:
  
  /// Default constructor
  HHOST_DEV CFMatSlice(T* data, size_t totNbRows, size_t totNbCols) : 
    EETMAT(MatExprT,CFMatSlice,T,N,M)(this), 
    m_totNbRows(totNbRows), m_totNbCols(totNbCols), m_data(data) {}
    
  /// Copy constructor from other slice
  HHOST_DEV CFMatSlice(const CFMatSlice<T,N,M>& in) : 
    EETMAT(MatExprT,CFMatSlice,T,N,M)(this), 
    m_totNbRows(in.m_totNbRows), m_totNbCols(in.m_totNbCols), m_data(in.m_data) {}
    
  /// Default destructor
  HHOST_DEV ~CFMatSlice() {}
  
  /// Overloading of assignment operator(s)
#define MATSLICE_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
    HHOST_DEV const CFMatSlice<T,N,M>& operator __op__ (MEETYPEV(EXPR) expr) \
    {EMATLOOP(i,0,N,j,0,M, (*this)(i COMMA j) __op__ expr.at(i COMMA j)); return *this;}
MATSLICE_EASSIGN_OP(=)
MATSLICE_EASSIGN_OP(+=)
MATSLICE_EASSIGN_OP(-=)
MATSLICE_EASSIGN_OP(*=)
MATSLICE_EASSIGN_OP(/=)
#undef MATSLICE_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define MATSLICE_EASSIGN_OP_CONST(__op__) \
  HHOST_DEV const CFMatSlice<T,N,M>& operator __op__ (T expr) {EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
MATSLICE_EASSIGN_OP_CONST(=)
MATSLICE_EASSIGN_OP_CONST(+=)
MATSLICE_EASSIGN_OP_CONST(-=)
MATSLICE_EASSIGN_OP_CONST(*=)
MATSLICE_EASSIGN_OP_CONST(/=)
#undef MATSLICE_EASSIGN_OP_CONST
  
  /// Overloading of the assignment operator "=" with CFVecSlice
  HHOST_DEV const CFMatSlice<T,N,M>& operator= (const CFMatSlice<T,N,M>& other)
  {
    cf_assert(&other != this);
    const size_t ns = size(); 
    for (size_t i = 0; i < ns; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// Overloading of the operator"()" for assignment.
  HHOST_DEV T& operator() (size_t i, size_t j) 
  {assert(i < N); assert(j < M); return m_data[i*m_totNbCols + j];}
 
  /// Overloading of the operator"()" (doesn't allow assignment)
  HHOST_DEV T operator() (size_t i, size_t j) const 
  {assert(i < N); assert(j < M); return m_data[i*m_totNbCols + j];}
  
  /// Overloading of the "[]" operator for assignment (writing).
  /// @param iElem index  
  HHOST_DEV T& operator[] (size_t iElem) 
  {assert(iElem < size()); return operator()(iElem/M, iElem%M);}
 
  /// Overloading of the "[]" operator for assignment (reading only).
  /// @param iElem index
  HHOST_DEV T operator[] (size_t iElem) const 
  {assert(iElem < size()); return operator()(iElem/M, iElem%M);}
  
  /// Accessor to individual entry
  /// @param iElem index
  HHOST_DEV T at (size_t iElem) const  
  {assert(iElem < size()); return operator()(iElem/M, iElem%M);}
  
  /// Accessor to individual entry
  /// @param i row index
  /// @param j column index    
  HHOST_DEV T at (size_t i, size_t j) const {return operator()(i,j);}
  
  /// return the array size 
  HHOST_DEV size_t size() const {return N*M;}
 
  /// return the number of rows 
  HHOST_DEV size_t nbRows() const {return N;}
  
  /// return the number of columns 
  HHOST_DEV size_t nbCols() const {return M;}
  
  /// copy content of another array
  template <typename ARRAY>
  HHOST_DEV void copyFrom(const ARRAY& in){EVECLOOP(i, 0, size(), m_data[i] = in[i]);}
  
  /// copy content of another array
  template <typename ARRAY>
  HHOST_DEV void copyTo(ARRAY& in) {EVECLOOP(i, 0, size(), in[i] = m_data[i]);}
    
  /// @return the raw data
  HHOST_DEV T* ptr() {return m_data;}
  
  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const CFMatSlice<T,N,M>& A);
    
  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>> 
    LTGT (std::istream& in, CFMatSlice<T,N,M>& A);
  
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
class CFMatSlice<T,0,0> : public EETMAT(MatExprT,CFMatSlice,T,0,0) {
public:
  
  // Constructor from preallocated memory
  HHOST_DEV CFMatSlice(T* data, size_t totNbRows, size_t totNbCols, size_t ns, size_t ms) : 
    EETMAT(MatExprT,CFMatSlice,T,0,0)(this),  
    m_totNbRows(totNbRows), m_totNbCols(totNbCols), m_nrows(ns), m_ncols(ms), m_data(data) {}
    
  /// Copy constructor from other slice
  HHOST_DEV CFMatSlice(const CFMatSlice<T,0,0>& in) : 
    EETMAT(MatExprT,CFMatSlice,T,0,0)(this),  
    m_totNbRows(in.m_totNbRows), m_totNbCols(in.m_totNbCols), m_nrows(in.m_nrows), m_ncols(in.m_ncols), m_data(in.m_data) {}
    
  /// Default destructor
  HHOST_DEV ~CFMatSlice() {}
  
   /// Overloading of assignment operator(s)
#define MATSLICE0_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
    HHOST_DEV const CFMatSlice<T,0,0>& operator __op__ (MEETYPEV(EXPR) expr) \
    {EMATLOOP(i,0,EGETSIZER(ENMAX(0,EXPR::SIZE1)),j,0,EGETSIZEC(ENMAX(0,EXPR::SIZE2)), (*this)(i COMMA j) __op__ expr.at(i COMMA j)); return *this;}
MATSLICE0_EASSIGN_OP(=)
MATSLICE0_EASSIGN_OP(+=)
MATSLICE0_EASSIGN_OP(-=)
MATSLICE0_EASSIGN_OP(*=)
MATSLICE0_EASSIGN_OP(/=)
#undef MATSLICE0_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define MATSLICE0_EASSIGN_OP_CONST(__op__) \
  HHOST_DEV const CFMatSlice<T,0,0>& operator __op__ (T expr) {EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
MATSLICE0_EASSIGN_OP_CONST(=)
MATSLICE0_EASSIGN_OP_CONST(+=)
MATSLICE0_EASSIGN_OP_CONST(-=)
MATSLICE0_EASSIGN_OP_CONST(*=)
MATSLICE0_EASSIGN_OP_CONST(/=)
#undef MATSLICE0_EASSIGN_OP_CONST

  /// Overloading of the assignment operator "=" with CFVecSlice
  HHOST_DEV const CFMatSlice<T,0,0>& operator= (const CFMatSlice<T,0,0>& other)
  {
    cf_assert(&other != this);
    const size_t ns = size(); 
    for (size_t i = 0; i < ns; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// Overloading of the operator"()" for assignment.
  HHOST_DEV T& operator() (size_t i, size_t j)
  {assert(i < m_nrows); assert(j < m_ncols); return m_data[i*m_totNbCols + j];}

  /// Overloading of the operator"()" (doesn't allow assignment)
  HHOST_DEV T operator() (size_t i, size_t j) const
  {assert(i < m_nrows); assert(j < m_ncols); return m_data[i*m_totNbCols + j];}
  
  /// Overloading of the "[]" operator for assignment (writing).
  /// @param iElem index  
  HHOST_DEV T& operator[] (size_t iElem) 
  {assert(iElem < size()); return operator()(iElem/m_ncols, iElem%m_ncols);}
 
  /// Overloading of the "[]" operator for assignment (reading only).
  /// @param iElem index
  HHOST_DEV T operator[] (size_t iElem) const 
  {assert(iElem < size()); return operator()(iElem/m_ncols, iElem%m_ncols);}
  
  /// Accessor to individual entry
  /// @param iElem index
  HHOST_DEV T at (size_t iElem) const  
  {assert(iElem < size()); return operator()(iElem/m_ncols, iElem%m_ncols);}
  
  /// Accessor to individual entry
  /// @param i row index
  /// @param j column index    
  HHOST_DEV T at (size_t i, size_t j) const {return operator()(i,j);}
  
  /// return the array size 
  HHOST_DEV size_t size() const {return m_nrows*m_ncols;}
 
  /// return the number of rows 
  HHOST_DEV size_t nbRows() const {return m_nrows;}
  
  /// return the number of columns 
  HHOST_DEV size_t nbCols() const {return m_ncols;}
  
  /// @return the raw data
  HHOST_DEV T* ptr() {return m_data;}
  
  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const CFMatSlice<T,0,0>& A);
    
  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>> 
    LTGT (std::istream& in, CFMatSlice<T,0,0>& A);
    
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

template <typename T, int N, int M>
std::ostream& operator<< (std::ostream& out, const CFMatSlice<T,N,M>& A)
{
  for (size_t i = 0; i < A.nbRows(); ++i) {
    for (size_t j = 0; j < A.nbCols(); ++j) {
      out << A(i,j) << " " ;
    }
    out << "\n";
  }
  return out;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N, int M>
std::istream& operator>> (std::istream& in, CFMatSlice<T,N,M>& A)
{
  for (size_t i = 0; i < A.nbRows(); ++i) {
    for (size_t j = 0; j < A.nbCols(); ++j) {
      in >> A(i,j);
    }
  }
  return in;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
