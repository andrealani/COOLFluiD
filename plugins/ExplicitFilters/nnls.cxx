// 	$Id: nnls.cc,v 1.1.2.1 2003/03/06 16:28:07 suvrit Exp $	
// File: nnls.cc
// Implements the Lawson-Hanson NNLS algorithm
// Copied over from nnls.c so i don't ahve copyright on this
//...somebody else has.....don't know if this should be under GPL or LGPL or
// whatever, but i'll put a lil note here anyways:

// Copyright (C) 2004 Lawson-Hanson
// Modifications to adapt to c++ by Suvrit Sra.

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#include "nnls.hh"
#include "MathTools/MathFunctions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////
 
void NonNegativeLeastSquares::solve(RealMatrix& A, RealVector& b, RealVector& x)
{
  // Convert variables for nhls function
  CFuint m=A.nbRows();
  CFuint n=A.nbCols();
  CFreal* a = new CFreal [m*n];
  CFreal* rhs = new CFreal [m];
  CFreal* solution = new CFreal[n];
  CFreal* w = new CFreal [n];
  CFreal* zz = new CFreal [n];
  CFuint* indx = new CFuint [n];
  for(CFuint i=0; i<m; ++i) {
    for (CFuint j=0; j<n; ++j) {
      a[i+j*m] = A(i,j);
    }
    rhs[i] = b[i];
  }
  CFuint mda = m;
  CFreal rnorm;
  CFuint mode;
  
  // Perform least squares
  nnls(a, mda, m, n, rhs, solution, &rnorm, w, zz, indx, &mode);
  cf_assert(mode==1);
  
  // Convert solution to COOLFluiD
  for (CFuint i=0; i<n; ++i) {
    x[i] = solution[i];
  }
  
  // Deallocations
  delete[] a;
  delete[] rhs;
  delete[] solution;
  delete[] w;
  delete[] zz;
  delete[] indx;
  
}

//////////////////////////////////////////////////////////////////////////////
 
CFuint NonNegativeLeastSquares::nnls(CFreal* a,  CFuint mda,  CFuint m,  CFuint n, CFreal* b, 
	 CFreal* x, CFreal* rnorm, CFreal* w, CFreal* zz, CFuint* index, 
	 CFuint* mode)
{
  /* System generated locals */
  CFuint a_dim1, a_offset, idx1, idx2;
  CFreal d1, d2;
 
 
  /* Local variables */
  static CFuint iter;
  static CFreal temp, wmax;
  static CFuint i__, j, l;
  static CFreal t, alpha, asave;
  static CFuint itmax, izmax, nsetp;
  static CFreal unorm, ztest, cc;
  CFreal dummy[2];
  static CFuint ii, jj, ip;
  static CFreal sm;
  static CFuint iz, jz;
  static CFreal up, ss;
  static CFuint rtnkey, iz1, iz2, npp1;
 
  /*     ------------------------------------------------------------------ 
   */
  /*     CFuint INDEX(N) */
  /*     CFreal precision A(MDA,N), B(M), W(N), X(N), ZZ(M) */
  /*     ------------------------------------------------------------------ 
   */
  /* Parameter adjustments */
  a_dim1 = mda;
  a_offset = a_dim1 + 1;
  a -= a_offset;
  --b;
  --x;
  --w;
  --zz;
  --index;
 
  /* Function Body */
  *mode = 1;
  if (m <= 0 || n <= 0) {
    *mode = 2;
    return 0;
  }
  iter = 0;
  itmax = n * 3;
 
  /*                    INITIALIZE THE ARRAYS INDEX() AND X(). */
 
  idx1 = n;
  for (i__ = 1; i__ <= idx1; ++i__) {
    x[i__] = 0.;
    /* L20: */
    index[i__] = i__;
  }
 
  iz2 = n;
  iz1 = 1;
  nsetp = 0;
  npp1 = 1;
  /*                             ******  MAIN LOOP BEGINS HERE  ****** */
 L30:
  /*                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION. 
   */
  /*                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED. */
 
  if (iz1 > iz2 || nsetp >= m) {
    goto L350;
  }
 
  /*         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W(). 
   */
 
  idx1 = iz2;
  for (iz = iz1; iz <= idx1; ++iz) {
    j = index[iz];
    sm = 0.;
    idx2 = m;
    for (l = npp1; l <= idx2; ++l) {
      /* L40: */
      sm += a[l + j * a_dim1] * b[l];
    }
    w[j] = sm;
    /* L50: */
  }
  /*                                   FIND LARGEST POSITIVE W(J). */
 L60:
  wmax = 0.;
  idx1 = iz2;
  for (iz = iz1; iz <= idx1; ++iz) {
    j = index[iz];
    if (w[j] > wmax) {
      wmax = w[j];
      izmax = iz;
    }
    /* L70: */
  }
 
  /*             IF WMAX .LE. 0. GO TO TERMINATION. */
  /*             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS. 
   */
 
  if (wmax <= 0.) {
    goto L350;
  }
  iz = izmax;
  j = index[iz];
 
  /*     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P. */
  /*     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID */
  /*     NEAR LINEAR DEPENDENCE. */
 
  asave = a[npp1 + j * a_dim1];
  idx1 = npp1 + 1;
  h12(c__1, &npp1, &idx1, m, &a[j * a_dim1 + 1], &c__1, &up, dummy, &
      c__1, &c__1, &c__0);
  unorm = 0.;
  if (nsetp != 0) {
    idx1 = nsetp;
    for (l = 1; l <= idx1; ++l) {
      /* L90: */
      /* Computing 2nd power */
      d1 = a[l + j * a_dim1];
      unorm += d1 * d1;
    }
  }
  unorm = std::sqrt(unorm);
  d2 = unorm + (d1 = a[npp1 + j * a_dim1], std::abs(d1)) * .01;
  if ((d2- unorm) > 0.) {
 
    /*        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE Z
	      Z */
    /*        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ). */
 
    idx1 = m;
    for (l = 1; l <= idx1; ++l) {
      /* L120: */
      zz[l] = b[l];
    }
    idx1 = npp1 + 1;
    h12(c__2, &npp1, &idx1, m, &a[j * a_dim1 + 1], &c__1, &up, (zz+1), &
	c__1, &c__1, &c__1);
    ztest = zz[npp1] / a[npp1 + j * a_dim1];
 
    /*                                     SEE IF ZTEST IS POSITIVE */
 
    if (ztest > 0.) {
      goto L140;
    }
  }
 
  /*     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P. */
  /*     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL */
  /*     COEFFS AGAIN. */
 
  a[npp1 + j * a_dim1] = asave;
  w[j] = 0.;
  goto L60;
 
  /*     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM */
  /*     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER */
  /*     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN */
  /*     COL J,  SET W(J)=0. */
 
 L140:
  idx1 = m;
  for (l = 1; l <= idx1; ++l) {
    /* L150: */
    b[l] = zz[l];
  }
 
  index[iz] = index[iz1];
  index[iz1] = j;
  ++iz1;
  nsetp = npp1;
  ++npp1;
 
  if (iz1 <= iz2) {
    idx1 = iz2;
    for (jz = iz1; jz <= idx1; ++jz) {
      jj = index[jz];
      h12(c__2, &nsetp, &npp1, m, 
	  &a[j * a_dim1 + 1], &c__1, &up, 
	  &a[jj * a_dim1 + 1], &c__1, &mda, &c__1);
      /* L160: */
    }
  }
 
  if (nsetp != m) {
    idx1 = m;
    for (l = npp1; l <= idx1; ++l) {
      /* L180: */
      // SS: CHECK THIS DAMAGE....
      a[l + j * a_dim1] = 0.;
    }
  }
 
  w[j] = 0.;
  /*                                SOLVE THE TRIANGULAR SYSTEM. */
  /*                                STORE THE SOLUTION TEMPORARILY IN ZZ(). 
   */
  rtnkey = 1;
  goto L400;
 L200:
 
  /*                       ******  SECONDARY LOOP BEGINS HERE ****** */
 
  /*                          ITERATION COUNTER. */
 
 L210:
  ++iter;
  if (iter > itmax) {
    *mode = 3;
    /* The following lines were replaced after the f2c translation */
    /* s_wsfe(&io___22); */
    /* do_fio(&c__1, " NNLS quitting on iteration count.", 34L); */
    /* e_wsfe(); */
    // CFlog(INFO, "\n NNLS quitting on iteration count.\n");
    goto L350;
  }
 
  /*                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE. */
  /*                                  IF NOT COMPUTE ALPHA. */
 
  alpha = 2.;
  idx1 = nsetp;
  for (ip = 1; ip <= idx1; ++ip) {
    l = index[ip];
    if (zz[ip] <= 0.) {
      t = -x[l] / (zz[ip] - x[l]);
      if (alpha > t) {
	alpha = t;
	jj = ip;
      }
    }
    /* L240: */
  }
 
  /*          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL */
  /*          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP. */
 
  if (alpha == 2.) {
    goto L330;
  }
 
  /*          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO */
  /*          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ. */
 
  idx1 = nsetp;
  for (ip = 1; ip <= idx1; ++ip) {
    l = index[ip];
    x[l] += alpha * (zz[ip] - x[l]);
    /* L250: */
  }
 
  /*        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I */
  /*        FROM SET P TO SET Z. */
 
  i__ = index[jj];
 L260:
  x[i__] = 0.;
 
  if (jj != nsetp) {
    ++jj;
    idx1 = nsetp;
    for (j = jj; j <= idx1; ++j) {
      ii = index[j];
      index[j - 1] = ii;
      g1(&a[j - 1 + ii * a_dim1], &a[j + ii * a_dim1], 
	 &cc, &ss, &a[j - 1 + ii * a_dim1]);
      // SS: CHECK THIS DAMAGE...
      a[j + ii * a_dim1] = 0.;
      idx2 = n;
      for (l = 1; l <= idx2; ++l) {
	if (l != ii) {
 
	  /*                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,
			     L)) */
 
	  temp = a[j - 1 + l * a_dim1];
	  // SS: CHECK THIS DAMAGE
	  a[j - 1 + l * a_dim1] = cc * temp + ss * a[j + l * a_dim1];
	  a[j + l * a_dim1] = -ss * temp + cc * a[j + l * a_dim1];
	}
	/* L270: */
      }
 
      /*                 Apply procedure G2 (CC,SS,B(J-1),B(J)) */
 
      temp = b[j - 1];
      b[j - 1] = cc * temp + ss * b[j];
      b[j] = -ss * temp + cc * b[j];
      /* L280: */
    }
  }
 
  npp1 = nsetp;
  --nsetp;
  --iz1;
  index[iz1] = i__;
 
  /*        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD 
   */
  /*        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED. */
  /*        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY */
  /*        THAT ARE NONPOSITIVE WILL BE SET TO ZERO */
  /*        AND MOVED FROM SET P TO SET Z. */
 
  idx1 = nsetp;
  for (jj = 1; jj <= idx1; ++jj) {
    i__ = index[jj];
    if (x[i__] <= 0.) {
      goto L260;
    }
    /* L300: */
  }
 
  /*         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK. */
 
  idx1 = m;
  for (i__ = 1; i__ <= idx1; ++i__) {
    /* L310: */
    zz[i__] = b[i__];
  }
  rtnkey = 2;
  goto L400;
 L320:
  goto L210;
  /*                      ******  END OF SECONDARY LOOP  ****** */
 
 L330:
  idx1 = nsetp;
  for (ip = 1; ip <= idx1; ++ip) {
    i__ = index[ip];
    /* L340: */
    x[i__] = zz[ip];
  }
  /*        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING. */
  goto L30;
 
  /*                        ******  END OF MAIN LOOP  ****** */
 
  /*                        COME TO HERE FOR TERMINATION. */
  /*                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR. */
 
 L350:
  sm = 0.;
  if (npp1 <= m) {
    idx1 = m;
    for (i__ = npp1; i__ <= idx1; ++i__) {
      /* L360: */
      /* Computing 2nd power */
      d1 = b[i__];
      sm += d1 * d1;
    }
  } else {
    idx1 = n;
    for (j = 1; j <= idx1; ++j) {
      /* L380: */
      w[j] = 0.;
    }
  }
  *rnorm = std::sqrt(sm);
  return 0;
 
  /*     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE */
  /*     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ(). */
 
 L400:
  idx1 = nsetp;
  for (l = 1; l <= idx1; ++l) {
    ip = nsetp + 1 - l;
    if (l != 1) {
      idx2 = ip;
      for (ii = 1; ii <= idx2; ++ii) {
	zz[ii] -= a[ii + jj * a_dim1] * zz[ip + 1];
	/* L410: */
      }
    }
    jj = index[ip];
    zz[ip] /= a[ip + jj * a_dim1];
    /* L430: */
  }
  switch ((CFuint)rtnkey) {
  case 1:  goto L200;
  case 2:  goto L320;
  }
 
  /* The next line was added after the f2c translation to keep
     compilers from complaining about a void return from a non-void
     function. */
  return 0;
 
} /* nnls_ */

//////////////////////////////////////////////////////////////////////////////

CFuint NonNegativeLeastSquares::g1(CFreal* a, CFreal* b, CFreal* cterm, CFreal* sterm, CFreal* sig)
{
  /* System generated locals */
  CFreal d;
 
  static CFreal xr, yr;
 
 
  if (std::abs(*a) > std::abs(*b)) {
    xr = *b / *a;
    /* Computing 2nd power */
    d = xr;
    yr = std::sqrt(d * d + 1.);
    d = 1. / yr;
    *cterm = MathTools::MathFunctions::changeSign(d, *a);
    *sterm = *cterm * xr;
    *sig = std::abs(*a) * yr;
    return 0;
  }
  if (*b != 0.) {
    xr = *a / *b;
    /* Computing 2nd power */
    d = xr;
    yr = std::sqrt(d * d + 1.);
    d = 1. / yr;
    *sterm = MathTools::MathFunctions::changeSign(d, *b);
    *cterm = *sterm * xr;
    *sig = std::abs(*b) * yr;
    return 0;
  }
  *sig = 0.;
  *cterm = 0.;
  *sterm = 1.;
  return 0;
} /* g1_ */
 
//////////////////////////////////////////////////////////////////////////////

/* See nnls.h for explanation */
CFuint NonNegativeLeastSquares::h12(CFuint mode, CFuint* lpivot, CFuint* l1, 
	CFuint m, CFreal* u, CFuint* iue, CFreal* up, CFreal* c__, 
	CFuint* ice, CFuint* icv, CFuint* ncv)
{
  /* System generated locals */
  CFuint u_dim1, u_offset, idx1, idx2;
  CFreal d, d2;
 
  /* Builtin functions */
  /* The following line was commented out after the f2c translation */
  /* CFreal std::sqrt(); */
 
  /* Local variables */
  static CFuint incr;
  static CFreal b;
  static CFuint i__, j;
  static CFreal clinv;
  static CFuint i2, i3, i4;
  static CFreal cl, sm;
 
  /*     ------------------------------------------------------------------ 
   */
  /*     CFreal precision U(IUE,M) */
  /*     ------------------------------------------------------------------ 
   */
  /* Parameter adjustments */
  u_dim1 = *iue;
  u_offset = u_dim1 + 1;
  u -= u_offset;
  --c__;
 
  /* Function Body */
  if (0 >= *lpivot || *lpivot >= *l1 || *l1 > m) {
    return 0;
  }
  cl = (d = u[*lpivot * u_dim1 + 1], std::abs(d));
  if (mode == 2) {
    goto L60;
  }
  /*                            ****** CONSTRUCT THE TRANSFORMATION. ****** 
   */
  idx1 = m;
  for (j = *l1; j <= idx1; ++j) {
    /* L10: */
    /* Computing MAX */
    d2 = (d = u[j * u_dim1 + 1], std::abs(d));
    cl = std::max(d2,cl);
  }
  if (cl <= 0.) {
    goto L130;
  } else {
    goto L20;
  }
 L20:
  clinv = 1. / cl;
  /* Computing 2nd power */
  d = u[*lpivot * u_dim1 + 1] * clinv;
  sm = d * d;
  idx1 = m;
  for (j = *l1; j <= idx1; ++j) {
    /* L30: */
    /* Computing 2nd power */
    d = u[j * u_dim1 + 1] * clinv;
    sm += d * d;
  }
  cl *= std::sqrt(sm);
  if (u[*lpivot * u_dim1 + 1] <= 0.) {
    goto L50;
  } else {
    goto L40;
  }
 L40:
  cl = -cl;
 L50:
  *up = u[*lpivot * u_dim1 + 1] - cl;
  u[*lpivot * u_dim1 + 1] = cl;
  goto L70;
  /*            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ****** 
   */
 
 L60:
  if (cl <= 0.) {
    goto L130;
  } else {
    goto L70;
  }
 L70:
  if (*ncv <= 0) {
    return 0;
  }
  b = *up * u[*lpivot * u_dim1 + 1];
  /*                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN. 
   */
 
  if (b >= 0.) {
    goto L130;
  } else {
    goto L80;
  }
 L80:
  b = 1. / b;
  i2 = 1 - *icv + *ice * (*lpivot - 1);
  incr = *ice * (*l1 - *lpivot);
  idx1 = *ncv;
  for (j = 1; j <= idx1; ++j) {
    i2 += *icv;
    i3 = i2 + incr;
    i4 = i3;
    sm = c__[i2] * *up;
    idx2 = m;
    for (i__ = *l1; i__ <= idx2; ++i__) {
      sm += c__[i3] * u[i__ * u_dim1 + 1];
      /* L90: */
      i3 += *ice;
    }
    if (sm != 0.) {
      goto L100;
    } else {
      goto L120;
    }
  L100:
    sm *= b;
    c__[i2] += sm * *up;
    idx2 = m;
    for (i__ = *l1; i__ <= idx2; ++i__) {
      c__[i4] += sm * u[i__ * u_dim1 + 1];
      /* L110: */
      i4 += *ice;
    }
  L120:
    ;
  }
 L130:
  return 0;
} /* h12 */
 
//////////////////////////////////////////////////////////////////////////////

  } // end namespace MathTools

} // end namespace COOLFluiD
