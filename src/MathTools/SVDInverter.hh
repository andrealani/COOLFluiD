// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_SVDInverter_hh
#define COOLFluiD_MathTools_SVDInverter_hh

//////////////////////////////////////////////////////////////////////////////

#include "RealVector.hh"
#include "MatrixInverter.hh"
#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a solver that uses the SVD decomposition,
 * based on Numerical Recipes v3.01
 * A = U*S*V'  -->  pinv(A) = V*inv(S)*U' 
 * U and V are orthogonal --> inv(U)=U' , inv(V)=V'
 * S is diagonal
 *
 * @author Willem Deconinck
 */
class MathTools_API SVDInverter : public MatrixInverter {
public:

  /**
   * Constructor
   *
   * @param   nbRows   number of rows
   * @param   nbCols   number of columns
   */
  SVDInverter(const CFuint& nbRows, const CFuint& nbCols);
  
  /**
   * Constructor
   * 
   * @param    a   The matrix to apply SVD to
   */
  SVDInverter(const RealMatrix& a);
  

  /**
   * Default destructor
   */
  ~SVDInverter();

  /** 
   * Invert the given matrix a and put the result in x
   * 
   * @param    a   The matrix to invert
   * @param    x   The inverted matrix
   */
  void invert(const RealMatrix& a, RealMatrix& x);
  
  /** 
   * Invert the given matrix a and put the result in x
   * 
   * @param    x   The inverted matrix
   */
  void invert(RealMatrix& x);

  /** 
   * Solve A x = b for a vector x using the pseudoinverse of A 
   * as obtained by SVD. If positive, thresh is the threshold 
   * value below which singular values are considered as zero. 
   * If thresh is negative, a default based on expected roundoﬀ 
   * error is used.
   * 
   * @param    a   The matrix of the system
   * @param    b   The right hand side of the system
   * @param    x   The solution of the system
   * @param    tresh   treshold value below which singular values are considered as zero
   */
	void solve(const RealMatrix& a, const RealVector& b, RealVector& x, CFreal thresh = -1.);
	
	/**
   * Solve a system A*x=b for x
   * with multiple right hand sides
   * 
   * @param    a   The matrix of the system
   * @param    b   The right hand sides of the system
   * @param    x   The solutions of the system
   * @param    tresh   treshold value below which singular values are considered as zero
   */
	void solve(const RealMatrix& a, const RealMatrix& b, RealMatrix& x, CFreal thresh = -1.);

  /** 
   * Solve A x = b for a vector x using the pseudoinverse of A 
   * as obtained by SVD. If positive, thresh is the threshold 
   * value below which singular values are considered as zero. 
   * If thresh is negative, a default based on expected roundoﬀ 
   * error is used.
   * 
   * @param    b   The right hand side of the system
   * @param    x   The solution of the system
   * @param    tresh   treshold value below which singular values are considered as zero
   */
   void solve(const RealVector& b, RealVector& x, CFreal thresh = -1.);
	
	/**
   * Solve a system A*x=b for x
   * with multiple right hand sides
   * 
   * @param    b   The right hand sides of the system
   * @param    x   The solutions of the system
   * @param    tresh   treshold value below which singular values are considered as zero
   */
	void solve(const RealMatrix& b, RealMatrix& x, CFreal thresh = -1.);

  /**
   * Get the matrix U where A=U*S*V'
   * @return   The matrix U
   */
  RealMatrix getU() { return _u;}

  /**
   * Get the matrix V where A=U*S*V'
   * @return   The matrix V
   */
  RealMatrix getV() { return _v;}
  
  /**
   * Get the matrix S where A=U*S*V'
   * @return   The singular values S (in vector)
   */
  RealVector getS() { return _s;}
  
  /**
   * @param    tresh   treshold value below which singular values are considered as zero
   * @return   The rank of the SVD problem
   */
  CFuint rank(CFreal thresh);
  
  /**
   * @param    tresh   treshold value below which singular values are considered as zero
   * @return   The nullity of the SVD problem
   */
  CFuint nullity(CFreal thresh);
  
  /**
   * @param    tresh   treshold value below which singular values are considered as zero
   * @return   The range of the SVD problem
   */
  RealMatrix range(CFreal thresh);
  
  /**
   * @param    tresh   treshold value below which singular values are considered as zero
   * @return   The nullspace of the SVD problem
   */
  RealMatrix nullspace(CFreal thresh);

  /**
   * The inverse of the condition number
   * @return   The inverse of the condition number
   */
  CFreal inv_condition() {
  	return (_s[0] <= 0. || _s[_cols-1] <= 0.) ? 0. : _s[_cols-1]/_s[0];
  }

  /**
   * Singular Value Decomposition
   */
  void decompose();
  
  /**
   * Reorder singular values in descending order
   */
  void reorder();

private: //helper functions

  CFreal pythag(const CFreal a, const CFreal b);
  
private: //data

  /// The number of rows and columns of A=U*S*V'
  CFuint                _rows, _cols;
  
  /// The matrices U and V  
  RealMatrix             _u, _v;

  /// The diagonal matrix S (in vector)
  RealVector             _s;
  
  /// Check to see if decomposed
  bool                   _decomposed;

}; // end of class SVDInverter

//////////////////////////////////////////////////////////////////////////////

class MathTools_API SVDException : public Common::Exception 
{
public:
  SVDException(const Common::CodeLocation& where, const std::string& what)
    : Common::Exception(where, what, "Singular Value Decomposition Exception") {}
    
  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  SVDException(const SVDException& e) throw() : Common::Exception(e) {}
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_SVDInverter_hh
