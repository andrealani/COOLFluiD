// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "JacobiEigenSolver.hh"
#include "MatrixEigenSolver.hh"
#include "MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

#define SMTX_MAX_ROTATIONS 20


//////////////////////////////////////////////////////////////////////////////
inline void JacobiEigenSolver::RotateMatrix(
    RealMatrix& a,
    CFreal& g, CFreal& h, CFreal& s, CFreal& c, CFreal& tau,
    const CFuint& i, const CFuint& j, const CFuint& k, const CFuint& l
    )
{
  g = a(i,j);
  h = a(k,l);
  a(i,j) = g - s*( h+g*tau );
  a(k,l) = h + s*( g-h*tau );
}
//////////////////////////////////////////////////////////////////////////////

void JacobiEigenSolver::eigenCalc( RealMatrix& a, RealMatrix& r, RealVector& lam)
{
  //CHECK( size == mNCols);

  cf_assert( a.nbRows() == a.nbCols() );

  CFuint i, k, size = a.nbRows();
  CFreal tresh, theta, tau, t, sm, s, h, g, c, tmp;

  RealVector b(size);
  RealVector z(size);
  RealVector w(size);


  // initialize
  //mtxL.InitDelta( size);
  //mtxD.InitDelta( size);


  lam = RealVector( 0.0, size);


  r = RealMatrix( size, size, 0.0);
  for ( CFuint ind=0; ind<size; ++ind)
    r( ind, ind) = 1.0;



  for ( CFuint ip=0; ip<size; ++ip)
  {
    b[ip] = w[ip] = a(ip,ip);
    z[ip] = 0.0;
  }


  // begin rotation sequence
  for ( i=0; i<SMTX_MAX_ROTATIONS; ++i)
  {
    sm = 0.0;
    for ( CFuint ip=0; ip<size-1; ++ip)
    {
      for ( CFuint iq=ip+1; iq<size; ++iq)
      {
        sm += ::std::abs( a(ip,iq) );
      }
    }

    if ( sm == 0.0)
      break;

    if (i < 3)                                // first 3 sweeps
      tresh = 0.2*sm/(size*size);
    else
      tresh = 0.0;


    for ( CFint ip=0; ip< static_cast<CFint>(size-1); ++ip)
    {
      for ( CFint iq=ip+1; iq<static_cast<CFint>(size); ++iq)
      {
        g = 100.0 * ::std::abs( a(ip,iq) );

        // after 4 sweeps
        if ( i > 3 && ( ::std::abs( w[ip])+g) == ::std::abs( w[ip]) && (::std::abs( w[iq])+g) == ::std::abs( w[iq]))
        {
          a(ip,iq) = 0.0;
        }
        else if ( ::std::abs( a(ip,iq) ) > tresh)
        {
          h = w[iq] - w[ip];
          if ( ( ::std::abs(h)+g) == ::std::abs(h))
          {
            t = ( a(ip,iq) ) / h;
          }
          else
          {
            theta = 0.5*h / ( a(ip,iq) );
            t = 1.0 / ( ::std::abs(theta) + ::sqrt( 1.0+theta*theta) );
            if ( theta < 0.0)
            {
              t = -t;
            }
          }
          c = 1.0 / ::sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*a(ip,iq);
          z[ip] -= h;
          z[iq] += h;
          w[ip] -= h;
          w[iq] += h;
          a(ip,iq) = 0.0;

          // ip already shifted left by 1 unit
          for ( CFint j = 0; j <= static_cast<CFint>(ip-1); ++j)
          {
            RotateMatrix( a, g, h, s, c, tau, j, ip, j, iq );
          }

          // ip and iq already shifted left by 1 unit
          for ( CFint j = ip+1; j <= static_cast<CFint>(iq-1); ++j)
          {
            RotateMatrix( a, g, h, s, c, tau, ip, j, j, iq );
          }

          // iq already shifted left by 1 unit
          for ( CFint j=iq+1; j<static_cast<CFint>(size); ++j)
          {
            RotateMatrix( a, g, h, s, c, tau, ip, j, iq, j );
          }

          for ( CFint j=0; j<static_cast<CFint>(size); ++j)
          {
            RotateMatrix( r, g, h, s, c, tau, j, ip, j, iq );
          }

        }
      }
    }

    for ( CFuint ip=0; ip<size; ++ip)
    {
      b[ip] += z[ip];
      w[ip] = b[ip];
      z[ip] = 0.0;
    }
  }


//   // this is NEVER called
//   if ( i >= SMTX_MAX_ROTATIONS )
//   {
//     std::cout << "SMatrix Decompose error: Error extracting eigenfunctions" << std::endl;
//   }


  // sort eigenfunctions                 these changes do not affect accuracy
  for ( CFuint j=0; j<size-1; ++j)                  // boundary incorrect
  {
    k = j;
    tmp = w[k];
    for ( CFuint i=j+1; i<size; ++i)                // boundary incorrect, shifted already
    {
      if ( w[i] > tmp)                   // why exchage if same?
      {
        k = i;
        tmp = w[k];
      }
    }
    if ( k != j)
    {
      w[k] = w[j];
      w[j] = tmp;
      for ( CFuint i=0; i<size; i++)
      {
        tmp = r(i,j);
        r(i,j) = r(i,k);
        r(i,k) = tmp;
      }
    }
  }


  for ( CFuint ip=0; ip<size; ++ip)
    lam[ip] = w[ip];

  return;


  // insure eigenvector consistency (i.e., Jacobi can compute vectors that
  // are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
  // reek havoc in hyperstreamline/other stuff. We will select the most
  // positive eigenvector.
  CFint  ceil_half_n = (static_cast<CFint>(size) >> 1) + (static_cast<CFint>(size) & 1);
  CFint numPos;

  for ( CFint j=0; j<static_cast<CFint>(size); j++)
  {
    for ( numPos=0, i=0; i<size; ++i)
    {
      if ( r(i,j) >= 0.0 )
      {
        ++numPos;
      }
    }
  //    if ( numPos < ceil(double(n)/double(2.0)) )
    if ( numPos < ceil_half_n)
    {
      for( CFint i=0; i<static_cast<CFint>(size); ++i)
      {
        r(i,j) *= -1.0;
      }
    }
  }

  for ( CFuint ip=0; ip<size; ++ip)
    lam[ip] = w[ip];

}


//////////////////////////////////////////////////////////////////////////////
/*
template< class ELEM_TYPE, MGSize MAX_SIZE>
inline void SMatrix< ELEM_TYPE, MAX_SIZE>::RotateMatrix(
    SMatrix< ELEM_TYPE, MAX_SIZE>& a,
    ELEM_TYPE& g, ELEM_TYPE& h, ELEM_TYPE& s, ELEM_TYPE& c, ELEM_TYPE& tau,
    const MGInt& i, const MGInt& j, const MGInt& k, const MGInt& l
    )
{
  g = a(i,j);
  h = a(k,l);
  a(i,j) = g - s*( h+g*tau );
  a(k,l) = h + s*( g-h*tau );
}


template< class ELEM_TYPE, MGSize MAX_SIZE>
inline void SMatrix< ELEM_TYPE, MAX_SIZE>::Decompose( SMatrix< ELEM_TYPE, MAX_SIZE>& mtxD, SMatrix< ELEM_TYPE, MAX_SIZE>& mtxL)
{
  CHECK( size == mNCols);


  MGSize i, k;
  ELEM_TYPE tresh, theta, tau, t, sm, s, h, g, c, tmp;

  SVector< ELEM_TYPE, MAX_SIZE> b(size);
  SVector< ELEM_TYPE, MAX_SIZE> z(size);
  SVector< ELEM_TYPE, MAX_SIZE> w(size);


  // initialize
  mtxL.InitDelta( size);
  mtxD.InitDelta( size);

  for ( MGSize ip=0; ip<size; ++ip)
  {
    b(ip) = w(ip) = (*this)(ip,ip);
    z(ip) = 0.0;
  }

  // begin rotation sequence
  for ( i=0; i<SMTX_MAX_ROTATIONS; ++i)
  {
    sm = 0.0;
    for ( MGSize ip=0; ip<size-1; ++ip)
    {
      for ( MGSize iq=ip+1; iq<size; ++iq)
      {
        sm += ::std::abs( (*this)(ip,iq) );
      }
    }

    if ( sm == 0.0)
      break;

    if (i < 3)                                // first 3 sweeps
      tresh = 0.2*sm/(size*size);
    else
      tresh = 0.0;


    for ( MGInt ip=0; ip<(MGInt)size-1; ++ip)
    {
      for ( MGInt iq=ip+1; iq<(MGInt)size; ++iq)
      {
        g = 100.0 * ::std::abs( (*this)(ip,iq) );

        // after 4 sweeps
        if ( i > 3 && ( ::std::abs(w(ip))+g) == ::std::abs(w(ip)) && (::std::abs(w(iq))+g) == ::std::abs(w(iq)))
        {
          (*this)(ip,iq) = 0.0;
        }
        else if ( ::std::abs( (*this)(ip,iq) ) > tresh)
        {
          h = w(iq) - w(ip);
          if ( ( ::std::abs(h)+g) == ::std::abs(h))
          {
            t = ( (*this)(ip,iq) ) / h;
          }
          else
          {
            theta = 0.5*h / ( (*this)(ip,iq) );
            t = 1.0 / ( ::std::abs(theta) + ::sqrt( 1.0+theta*theta) );
            if ( theta < 0.0)
            {
              t = -t;
            }
          }
          c = 1.0 / ::sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*(*this)(ip,iq);
          z(ip) -= h;
          z(iq) += h;
          w(ip) -= h;
          w(iq) += h;
          (*this)(ip,iq) = 0.0;

          // ip already shifted left by 1 unit
          for ( MGInt j = 0; j <= (MGInt)ip-1; ++j)
          {
            RotateMatrix( *this, g, h, s, c, tau, j, ip, j, iq );
          }

          // ip and iq already shifted left by 1 unit
          for ( MGInt j = ip+1; j <= (MGInt)iq-1; ++j)
          {
            RotateMatrix( *this, g, h, s, c, tau, ip, j, j, iq );
          }

          // iq already shifted left by 1 unit
          for ( MGInt j=iq+1; j<(MGInt)size; ++j)
          {
            RotateMatrix( *this, g, h, s, c, tau, ip, j, iq, j );
          }

          for ( MGInt j=0; j<(MGInt)size; ++j)
          {
            RotateMatrix( mtxL, g, h, s, c, tau, j, ip, j, iq );
          }

        }
      }

    }

    for ( MGSize ip=0; ip<size; ++ip)
    {
      b(ip) += z(ip);
      w(ip) = b(ip);
      z(ip) = 0.0;
    }
  }


  // this is NEVER called
  if ( i >= SMTX_MAX_ROTATIONS )
  {
    THROW_SMATRIX( "SMatrix Decompose error: Error extracting eigenfunctions");
  }


  // sort eigenfunctions                 these changes do not affect accuracy
  for ( MGSize j=0; j<size-1; ++j)                  // boundary incorrect
  {
    k = j;
    tmp = w(k);
    for ( MGSize i=j+1; i<size; ++i)                // boundary incorrect, shifted already
    {
      if (w(i) >= tmp)                   // why exchage if same?
      {
        k = i;
        tmp = w(k);
      }
    }
    if ( k != j)
    {
      w(k) = w(j);
      w(j) = tmp;
      for ( MGSize i=0; i<size; i++)
      {
        tmp = mtxL(i,j);
        mtxL(i,j) = mtxL(i,k);
        mtxL(i,k) = tmp;
      }
    }
  }


  for ( MGSize ip=0; ip<size; ++ip)
    mtxD(ip,ip) = w(ip);
  return;


  // insure eigenvector consistency (i.e., Jacobi can compute vectors that
  // are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
  // reek havoc in hyperstreamline/other stuff. We will select the most
  // positive eigenvector.
  MGInt  ceil_half_n = ((MGInt)size >> 1) + ((MGInt)size & 1);
  MGInt numPos;

  for ( MGInt j=0; j<(MGInt)size; j++)
  {
    for ( numPos=0, i=0; i<size; ++i)
    {
      if ( mtxL(i,j) >= 0.0 )
      {
        ++numPos;
      }
    }
  //    if ( numPos < ceil(double(n)/double(2.0)) )
    if ( numPos < ceil_half_n)
    {
      for( MGInt i=0; i<(MGInt)size; ++i)
      {
        mtxL(i,j) *= -1.0;
      }
    }
  }

  for ( MGSize ip=0; ip<size; ++ip)
    mtxD(ip,ip) = w(ip);

}
*/

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
