// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LUInverter.hh"
#include "MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

LUInverter::LUInverter(const CFuint& n) :
  MatrixInverter(),
  _n(n),
  _indx(n),
  _tempcol(n),
  _vv(n),
  _a(n,n,0.0)
{
}

//////////////////////////////////////////////////////////////////////////////

LUInverter::~LUInverter()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUInverter::invert(const RealMatrix& a, RealMatrix& x)
{
  _a = a;
  factLU();

  for (CFuint j = 0;  j < _n; ++j) {
    _tempcol = 0.0;
    _tempcol[j] = 1.0;
    solveForwBack();
    for (CFuint i = 0; i < _n; ++i) {
      x(i,j) = _tempcol[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LUInverter::factLU()
{
  CFuint imax = 0;
  CFreal big = 0.0;
  CFreal dum = 0.0;
  CFreal sum = 0.0;
  CFreal temp = 0.0;
  CFreal d = 1.0;  // no row permutation yet

  for (CFuint i = 0; i < _n; ++i) {
    big = 0.0;
    for (CFuint j = 0; j < _n; ++j) {
      if ((temp = std::abs(_a(i,j))) > big) {
	big = temp;
      }
    }
    cf_assert(std::abs(big) > MathTools::MathConsts::CFrealEps());
    _vv[i] = 1.0/big;
  }
  for (CFuint j = 0; j < _n; ++j) {
    for (CFuint i = 0; i < j; ++i) {
      sum = _a(i,j);
      for (CFuint k = 0; k < i; ++k) {
	sum -= _a(i,k)*_a(k,j);
      }
      _a(i,j) = sum;
    }
    big = 0.0;

    for (CFuint i = j; i < _n; ++i) {
      sum = _a(i,j);
      for (CFuint k = 0; k < j; ++k) {
	sum -= _a(i,k)*_a(k,j);
      }
      _a(i,j) = sum;
      if ((dum = _vv[i]*std::abs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (CFuint k = 0; k < _n; ++k) {
	dum = _a(imax,k);
	_a(imax,k) = _a(j,k);
	_a(j,k) = dum;
      }
      d = -d;
      _vv[imax] = _vv[j];
    }
    _indx[j] = imax;
    cf_assert(MathChecks::isNotZero(_a(j,j)));
    if (j != _n) {
      dum = 1.0 / _a(j,j);
      for (CFuint i = j+1; i < _n; ++i) {
	_a(i,j) *= dum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LUInverter::solveForwBack()
{
  CFint ip = 0;
  CFint ii = -1;
  CFreal sum = 0.0;

  for (CFuint i = 0; i < _n; ++i) {
    ip = _indx[i];
    sum = _tempcol[ip];
    _tempcol[ip] = _tempcol[i];
    if (ii >= 0) {
      for (CFuint j = ii; j <= i-1; ++j) {
	sum -= _a(i,j)*_tempcol[j];
      }
    }
    else if (sum) {
      ii = i;
    }
    _tempcol[i] = sum;
  }
  for (CFint i = _n-1; i >= 0; --i) {
    sum = _tempcol[i];
    for (CFuint j = i+1; j < _n; ++j) {
      sum -= _a(i,j)*_tempcol[j];
    }
    _tempcol[i] = sum/_a(i,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
