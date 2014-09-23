// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_LUInverterT_hh
#define COOLFluiD_MathTools_LUInverterT_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// This class inverts a matrix using LU decomposition
/// @author Tiago Quintino
template < unsigned int SIZE >
struct LUInverterT {
public:

  /// Constructor
  LUInverterT ();

  /// Invert the given matrix a and put the result in x
  void invert (const RealMatrix& a, RealMatrix& x);

  /// Factorize the matrix
  void factorizeLU();

  /// Solve with forward and backward substitution
  void solveForwBack();

private: // data
  /// storage of indexes
  std::valarray<CFuint> m_indx;
  /// temporary vector
  RealVector              m_tmp_col;
  /// temporary vector
  RealVector              m_vv;
  /// temporary copy of input matrix
  RealMatrix              m_a;

}; // end of class LUInverter

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/LUInverterT.ci"

//////////////////////////////////////////////////////////////////////////////
#endif // COOLFluiD_MathTools_LUInverterT_hh
