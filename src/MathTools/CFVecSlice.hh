// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_CFVecSlice_hh
#define COOLFluiD_MathTools_CFVecSlice_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/ExprT.hh"
#include "Common/Compatibility.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
  
  template  <typename T, int N> class CFVecSlice;
  
  template <typename T, int N> std::ostream& operator<< 
    (std::ostream& out, const CFVecSlice<T,N>& v);
  
  template <typename T, int N> std::istream& operator>> 
    (std::istream& in, CFVecSlice<T,N>& v);
  
//////////////////////////////////////////////////////////////////////////////
  
/// Definition of a class CFVecSlice that implements an expression template technique
/// with statically defined size
/// @author Andrea Lani
template <typename T, int N = 0>
class CFVecSlice : public EETVEC(ExprT,CFVecSlice,T,N) {
public:
  
  /// Default constructor
  HHOST_DEV CFVecSlice(T* data) : EETVEC(ExprT,CFVecSlice,T,N)(this), m_data(data) {}
    
  /// Copy constructor from other slice
  HHOST_DEV CFVecSlice(const CFVecSlice<T,N>& in) : EETVEC(ExprT,CFVecSlice,T,N)(this), m_data(in.m_data) {}
      
  /// Default destructor
  HHOST_DEV ~CFVecSlice() {}
  
  /// Overloading of assignment operator(s)
#define VECSLICE_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
  HHOST_DEV const CFVecSlice<T,N>& operator __op__ (EETYPEV(EXPR) expr) \
    {EVECLOOP(i, 0, N, m_data[i] __op__ expr.at(i)); return *this;}
VECSLICE_EASSIGN_OP(=)
VECSLICE_EASSIGN_OP(+=)
VECSLICE_EASSIGN_OP(-=)
VECSLICE_EASSIGN_OP(*=)
VECSLICE_EASSIGN_OP(/=)
#undef VECSLICE_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define VECSLICE_EASSIGN_OP_CONST(__op__) \
  HHOST_DEV const CFVecSlice<T,N>& operator __op__ (T expr) \
  {EVECLOOP(i, 0, N, m_data[i] __op__ expr); return *this;}
VECSLICE_EASSIGN_OP_CONST(=)
VECSLICE_EASSIGN_OP_CONST(+=)
VECSLICE_EASSIGN_OP_CONST(-=)
VECSLICE_EASSIGN_OP_CONST(*=)
VECSLICE_EASSIGN_OP_CONST(/=)
#undef VECSLICE_EASSIGN_OP_CONST

  /// Overloading of the assignment operator "=" with CFVecSlice
  HHOST_DEV const CFVecSlice<T,N>& operator= (const CFVecSlice<T,N>& other)
  {
    cf_assert(&other != this);
    for (size_t i = 0; i < N; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// return the array size 
  HHOST_DEV size_t size() const {return N;}

  /// copy content of another array
  template <typename ARRAY>
  HHOST_DEV void copyFrom(const ARRAY& in){EVECLOOP(i, 0, N, m_data[i] = in[i]);}
  
  /// copy content of another array
  template <typename ARRAY>
  HHOST_DEV void copyTo(ARRAY& in) {EVECLOOP(i, 0, N, in[i] = m_data[i]);}
    
  /// @return the raw data
  HHOST_DEV T* ptr() {return m_data;}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const CFVecSlice<T,N>& v);
    
  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>> 
    LTGT (std::istream& in, CFVecSlice<T,N>& v);
  
private:
  
  /// array data
  T* m_data;
};
  
//////////////////////////////////////////////////////////////////////////////

/// Definition of a class CFVecSlice that implements an expression template technique
/// with dynamical size
/// @author Andrea Lani
template <typename T>
class CFVecSlice<T,0> : public EETVEC(ExprT,CFVecSlice,T,0) {
public:
  
  // Constructor from preallocated memory
  HHOST_DEV CFVecSlice(T* data, size_t ns) : 
    EETVEC(ExprT,CFVecSlice,T,0)(this), m_size(ns), m_data(data) {}
    
  /// Copy constructor from other slice
  HHOST_DEV CFVecSlice(const CFVecSlice<T,0>& in) : 
    EETVEC(ExprT,CFVecSlice,T,0)(this), m_size(in.size()), m_data(in.m_data) {}
    
  /// Default destructor
  HHOST_DEV ~CFVecSlice() {}
  
   /// Overloading of assignment operator(s)
#define VECSLICE0_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
  HHOST_DEV const CFVecSlice<T,0>& operator __op__ (EETYPEV(EXPR) expr) \
    {EVECLOOP(i, 0, EGETSIZE1(ENMAX(0,EXPR::SIZE1)), m_data[i] __op__ expr.at(i)); return *this;}
VECSLICE0_EASSIGN_OP(=)
VECSLICE0_EASSIGN_OP(+=)
VECSLICE0_EASSIGN_OP(-=)
VECSLICE0_EASSIGN_OP(*=)
VECSLICE0_EASSIGN_OP(/=)
#undef VECSLICE0_EASSIGN_OP

/// Overloading of assignment operator(s) with constants
#define VECSLICE0_EASSIGN_OP_CONST(__op__) \
  HHOST_DEV const CFVecSlice<T,0>& operator __op__ (T expr) \
  {EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
VECSLICE0_EASSIGN_OP_CONST(=)
VECSLICE0_EASSIGN_OP_CONST(+=)
VECSLICE0_EASSIGN_OP_CONST(-=)
VECSLICE0_EASSIGN_OP_CONST(*=)
VECSLICE0_EASSIGN_OP_CONST(/=)
#undef VECSLICE0_EASSIGN_OP_CONST
 
  /// Overloading of the assignment operator "=" with CFVecSlice
  HHOST_DEV const CFVecSlice<T,0>& operator= (const CFVecSlice<T,0>& other)
  {
    cf_assert(&other != this);
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
 
  /// return the array size 
  HHOST_DEV size_t size() const {return m_size;}
 
  /// @return the raw data
  HHOST_DEV T* ptr() {return m_data;}
 
  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const CFVecSlice<T,0>& v);
    
  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>> 
    LTGT (std::istream& in, CFVecSlice<T,0>& v);
  
private:
  
  /// array size
  size_t m_size;
  
  /// array data
  T* m_data;
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
std::ostream& operator<< (std::ostream& out, const CFVecSlice<T,N>& v)
{
  const size_t size = v.size();
  for (size_t i = 0; i < size; ++i)
    out << v.m_data[i] << " " ;
  return out;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
std::istream& operator>> (std::istream& in, CFVecSlice<T,N>& v)
{
  const size_t size = v.size();
  for (size_t i = 0; i < size; ++i)
    in >> v.m_data[i];
  return in;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
