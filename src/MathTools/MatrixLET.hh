// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MatrixLET_hh
#define COOLFluiD_MathTools_MatrixLET_hh

//////////////////////////////////////////////////////////////////////////////

#include <ostream>

#include "MathTools/LExprOp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
    
    template <class T, class TAG, CFuint N> 
    class MatrixLET;
    
    template <class T, class TAG, CFuint N> 
    std::ostream& operator<< (std::ostream& out, const MatrixLET<T,TAG,N>& v);
    
    template <class T, class TAG, CFuint N>
    std::istream& operator>> (std::istream& in, MatrixLET<T,TAG,N>& v);
    
    template <class T, class TAG> 
    std::ostream& operator<< (std::ostream& out, const MatrixLET<T,TAG,0>& v);
    
    template <class T, class TAG>
    std::istream& operator>> (std::istream& in, MatrixLET<T,TAG,0>& v);
    
//////////////////////////////////////////////////////////////////////////////

/**
 * Definition of a high-performance matrix class for numerical applications 
 * that uses an improved version of the fast expression template concept.
 * This techique allows to retain a user-friendly Matlab-like syntax at the cost
 * of two extra template parameters: 
 * 1- a Tag class which ensures unicity of the vector instantiation within 
 *    the context of a class) 
 * 2- an integer that distinguish a C style allocated storage with fixed size from 
 *    a dynamical allocated one.
 *
 * @author Andrea Lani
 *
 */
template <class T, class TAG, CFuint N = 0>
class MatrixLET : public LExpr<MatrixLET<T,TAG,N>,T,N> {
public:
  /**
   * Default Constructor.
   */
  MatrixLET();
  
  /**
   * Destructor.
   */
  ~MatrixLET();
  
   /**
   * Overloading of the stream operator "<<" for the output
   * No "\n"ine introduced.
   *
   * @param out missing documentation
   * @param v missing documentation
   * @return missing documentation
   */
  friend std::ostream& operator<< LTGT (std::ostream& out, 
					const MatrixLET<T,TAG,N>& v);
  
  /**
   * Overloading of the stream operator ">>" for the input
   *
   * @param in missing documentation
   * @param v missing documentation
   * @return missing documentation
   */
  friend std::istream& operator>> LTGT (std::istream& in, 
					 MatrixLET<T,TAG,N>& v);
  
  /**
   * Accessor used by the expression template engine
   */
  static inline T at (CFuint i) 
  {
    return _v[i];
  }
  
  /**
   * Overloading of the "[]" operator for assignment (writing).
   * @param i index
   */
  T& operator[] (CFuint i) 
  {
    cf_assert(i < size());
    return _v[i];
  }
  
  /**
   * Overloading of the "[]" operator for assignment (writing).
   * @param i index
   */
  const T& operator[] (CFuint i) const
  {
    cf_assert(i < size());
    return _v[i];
  } 
  
  /**
   * Overloading of the "()" operator for assignment (writing).
   * @param i index
   * @param j index
   */
  T& operator() (CFuint i, CFuint j) 
  {
    cf_assert(i*N + j < size());
    return _v[i*N + j];
  }
  
  /**
   * Overloading of the "()" operator for assignment (reading).
   * @param i index
   * @param j index
   */
  const T& operator() (CFuint i, CFuint j) const 
  {
    cf_assert(i*N + j < size());
    return _v[i*N + j];
  } 
  
  /**
   * @return the size (the number of elements) of the MatrixLET
   */
  CFuint size() const  
  {
    return N*N;
  }
  
  /**
   * @return the number of rows in the MatrixLET
   */
  static inline CFuint nbRows() {return N;}
  
  /**
   * @return the number of columns in the MatrixLET
   */
  static inline CFuint nbCols() {return N;}
  
  /**
   * Overloading for operator= taking an expression as argument
   */ 
#define MATLET_EQ_OP(__op__)				\
template <class EXPR>							\
const MatrixLET& operator __op__ (const LExpr<EXPR,T,EXPR::SIZE>& expr)	\
{									\
  const CFuint nmax = N*N;						\
  for (CFuint i = 0; i < nmax; ++i) {					\
    _v[i] __op__ EXPR::at(i);						\
  }									\
  return *this;								\
}
  
  MATLET_EQ_OP(=)
  MATLET_EQ_OP(+=)
  MATLET_EQ_OP(-=)
    
#undef MATLET_EQ_OP

private:
  /// static array
  static T _v[N*N];
};

//////////////////////////////////////////////////////////////////////////////

/**
 * Partial specialization (supporting dynamical memory allocation) of the 
 * high-performance matrix class for numerical applications that uses 
 * an improved version of the fast expression template concept.
 * This techique allows to retain a user-friendly Matlab-like syntax at 
 * the cost of two extra template parameters: 
 * 1- a Tag class which ensures unicity of the vector instantiation within 
 *    the context of a class) 
 * 2- an integer that distinguish a C style allocated storage with fixed 
 *    size from a dynamical allocated one (here set == 0).
 *
 * @author Andrea Lani
 *
 */
template <class T, class TAG>
class MatrixLET<T,TAG,0> : public LExpr<MatrixLET<T,TAG,0>,T,0> {
public:
  
  /**
   * Default Constructor.
   */
  MatrixLET(CFuint nr, CFuint nc);
  
  /**
   * Destructor.
   */
  ~MatrixLET();
  
  /**
   * Overloading of operator= in order to convert a given vector of type OTHER
   * to a VectorLET
   */
  template <class OTHER>
  void setPtr(OTHER& array)
  {
    _v = &array[0];
  }

  /**
   * This function must be called ater having used the operator=(CFVector)
   * in order to reset the pointer to default behaviour
   */
  void release()
  {
    _nr = 0;
    _nc = 0;
    _v = CFNULL;
  }
  
   /**
   * Overloading of the stream operator "<<" for the output
   * No "\n"ine introduced.
   *
   * @param out missing documentation
   * @param v missing documentation
   * @return missing documentation
   */
  friend std::ostream& operator<< LTGT (std::ostream& out, 
					const MatrixLET<T,TAG,0>& v);
  
  /**
   * Overloading of the stream operator ">>" for the input
   *
   * @param in missing documentation
   * @param v missing documentation
   * @return missing documentation
   */
  friend std::istream& operator>> LTGT (std::istream& in, 
					 MatrixLET<T,TAG,0>& v);
  
  /**
   * Accessor used by the expression template engine
   */
  static inline T at (CFuint i) 
  {
    return _v[i];
  }
  
  /**
   * Overloading of the "[]" operator for assignment (writing).
   * @param i index
   */
  T& operator[] (CFuint i) 
  {
    cf_assert(i < size());
    return _v[i];
  }
  
  /**
   * Overloading of the "[]" operator for assignment (writing).
   * @param i index
   */
  const T& operator[] (CFuint i) const
  {
    cf_assert(i < size());
    return _v[i];
  } 
  
  /**
   * Overloading of the "()" operator for assignment (writing).
   * @param i index
   * @param j index
   */
  T& operator() (CFuint i, CFuint j) 
  {
    cf_assert(i*_nr + j < size());
    return _v[i*_nr + j];
  }
  
  /**
   * Overloading of the "()" operator for assignment (reading).
   * @param i index
   * @param j index
   */
  const T& operator() (CFuint i, CFuint j) const 
  {
    cf_assert(i*_nr + j < size());
    return _v[i*_nr + j];
  } 
  
  /**
   * @return the size (the number of elements) of the MatrixLET
   */
  CFuint size() const  
  {
    return _nr*_nc;
  }
  
  /**
   * @return the number of rows in the MatrixLET
   */
  static inline CFuint nbRows() {return _nr;}
  
  /**
   * @return the number of columns in the MatrixLET
   */
  static inline CFuint nbCols() {return _nc;}
  
  /**
   * Overloading for operator= taking an expression as argument
   */ 
#define MATLET_EQ_OP(__op__)				\
template <class EXPR>							\
const MatrixLET& operator __op__ (const LExpr<EXPR,T,EXPR::SIZE>& expr)	\
{									\
   const CFuint nmax = _nc*_nr;				\
    for (CFuint i = 0; i < nmax; ++i) {					\
    _v[i] __op__ EXPR::at(i);						\
  }									\
  return *this;								\
}
  
  MATLET_EQ_OP(=)
  MATLET_EQ_OP(+=)
  MATLET_EQ_OP(-=)
    
#undef MATLET_EQ_OP
    
private:
  /// number of rows
  static CFuint _nr;
  
  /// number of columns
  static CFuint _nc;
  
  /// static pointer to raw data
  static T*     _v;
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MatrixLET.ci"

//////////////////////////////////////////////////////////////////////////////

#endif
