// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "MathTools/IntersectSolver.hh"
#include "MathTools/MatrixIntersect.hh"
#include "MathTools/MatrixInverter.hh"
#include "MathTools/MatrixEigenSolver.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {



//////////////////////////////////////////////////////////////////////////////

void IntersectSolver::intersectCalc( const RealMatrix& a, const RealMatrix& b, RealMatrix& res )
{
  //const CFuint size = 3;
  const CFuint size = a.nbRows();

  cf_assert( a.nbRows() == a.nbCols() );
  cf_assert( b.nbRows() == b.nbCols() );
  cf_assert( size == a.nbRows() );
  cf_assert( size == b.nbRows() );

  MatrixEigenSolver *pEigenSol = MatrixEigenSolver::create( size, true);
  MatrixInverter *pInverter = MatrixInverter::create( size, false);


  RealMatrix mtx_a(a);
  RealMatrix mtxL_a( size, size, 0.0), mtxLI_a( size, size, 0.0);
  RealMatrix mtxND_a( size, size, 0.0);
  RealVector vecD_a( 0.0, size), vecND_a( 0.0, size);

  RealMatrix mtx_b(b);
  RealMatrix mtxL_b( size, size, 0.0), mtxLI_b( size, size, 0.0);
  RealMatrix mtxND_b( size, size, 0.0);
  RealVector vecD_b( 0.0, size), vecND_b( 0.0, size);


  pEigenSol->eigenCalc( mtx_a, mtxL_a, vecD_a);
  pInverter->invert( mtxL_a, mtxLI_a);

  pEigenSol->eigenCalc( mtx_b, mtxL_b, vecD_b);
  pInverter->invert( mtxL_b, mtxLI_b);


  mtxND_a = mtxLI_a * (vecD_b * mtxL_a);
  mtxND_b = mtxLI_b * (vecD_a * mtxL_b);

  for ( CFuint i=0; i<size; ++i)
  {
//    vecND_a[i] = max( mtxND_a(i,i), vecD_a[i] );
//    vecND_b[i] = max( mtxND_b(i,i), vecD_b[i] );
    vecND_a[i] = mtxND_a(i,i) > vecD_a[i]  ? mtxND_a(i,i) : vecD_a[i];
    vecND_b[i] = mtxND_b(i,i) > vecD_b[i]  ? mtxND_b(i,i) : vecD_b[i];
  }

  mtx_a = mtxL_a * (vecND_a * mtxLI_a);
  mtx_b = mtxL_b * (vecND_b * mtxLI_b);

  res = 0.5 * ( mtx_a + mtx_b );


  deletePtr( pInverter);
  deletePtr( pEigenSol);
}


//////////////////////////////////////////////////////////////////////////////


  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
