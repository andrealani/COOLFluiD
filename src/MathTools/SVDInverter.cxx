// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SVDInverter.hh"
#include "MathChecks.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

SVDInverter::SVDInverter(const CFuint& nbRows, const CFuint& nbCols) :
  MatrixInverter(),
  _rows(nbRows),
  _cols(nbCols),
  _u(nbRows,nbCols),
  _v(nbCols,nbCols),
  _s(nbCols),
  _decomposed(false)
{
}

//////////////////////////////////////////////////////////////////////////////

SVDInverter::SVDInverter(const RealMatrix& a) :
  MatrixInverter(),
  _rows(a.nbRows()),
  _cols(a.nbCols()),
  _u(a),
  _v(_cols,_cols),
  _s(_cols),
  _decomposed(false)
{
  try {
    decompose();
  }
  catch(SVDException& e) {
     CFLog(VERBOSE, e.what() << "\n");
     throw;
  }
  reorder();
}

//////////////////////////////////////////////////////////////////////////////

SVDInverter::~SVDInverter()
{
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::invert(const RealMatrix& a, RealMatrix& x)
{
  if (!_decomposed) {
    _u = a;
    decompose();
    reorder();
	}

  invert(x);
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::invert(RealMatrix& x)
{
  // cf_assert(_decomposed == true);
  
  CFreal tsh = 0.5*std::sqrt(_rows+_cols+1.)*_s[0]*MathConsts::CFrealEps();
  
  RealVector invS(_cols);  
  for (CFuint j=0; j<_cols; ++j) {
    if (_s[j] > tsh) {
		  invS[j]=1./_s[j];
	  } else {
      invS[j]=0.0;
	  }
  }
    
  RealMatrix transU(_cols,_rows);
  _u.transpose(transU);
  
  if (x.nbRows()!=_cols || x.nbCols()!=_rows){
    x.resize(_cols,_rows);
  }
  // Moore-Penrose Inverse  pinv(A) = V*inv(S)*U' 
  x = _v * (invS * transU);
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::solve(const RealMatrix& a, const RealVector& b, RealVector& x, CFreal thresh)
{  
  if (!_decomposed) {
    _u = a;
    decompose();
    reorder();
	}
  solve(b,x,thresh);
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::solve(const RealMatrix& a, const RealMatrix& b, RealMatrix& x, CFreal thresh)
{
	CFuint i,j,m=b.nbCols();
	if (b.nbRows() != _cols || x.nbRows() != _cols || b.nbCols() != x.nbCols())
		throw SVDException(FromHere(),"solve bad sizes");
	RealVector xx(_cols);
	for (j=0;j<m;j++) {
		for (i=0;i<_cols;i++) xx[i] = b(i,j);
		solve(a,xx,xx,thresh);
		for (i=0;i<_cols;i++) x(i,j) = xx[i];
	}
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::solve(const RealVector& b, RealVector& x, CFreal thresh) {
  
	CFuint i,j,jj;
	CFreal s;
	if (b.size() != _rows || x.size() != _cols) throw SVDException(FromHere(),"solve bad sizes");
	RealVector tmp(_cols);
	CFreal tsh = (thresh >= 0. ? thresh : 0.5*std::sqrt(_rows+_cols+1.)*_s[0]*MathConsts::CFrealEps());
	for (j=0;j<_cols;j++) {
		s=0.0;
		if (_s[j] > tsh) {
			for (i=0;i<_rows;i++) s += _u(i,j)*b[i];
			s /= _s[j];
		}
		tmp[j]=s;
	}
	for (j=0;j<_cols;j++) {
		s=0.0;
		for (jj=0;jj<_cols;jj++) s += _v(j,jj)*tmp[jj];
		x[j]=s;
	}
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::solve(const RealMatrix& b, RealMatrix& x, CFreal thresh)
{
	CFuint i,j,m=b.nbCols();
	if (b.nbRows() != _cols || x.nbRows() != _cols || b.nbCols() != x.nbCols())
		throw SVDException(FromHere(),"solve bad sizes");
	RealVector xx(_cols);
	for (j=0;j<m;j++) {
		for (i=0;i<_cols;i++) xx[i] = b(i,j);
		solve(xx,xx,thresh);
		for (i=0;i<_cols;i++) x(i,j) = xx[i];
	}
}

//////////////////////////////////////////////////////////////////////////////

CFuint SVDInverter::rank(CFreal thresh = -1.) {
	CFuint j,nr=0;
	CFreal tsh = (thresh >= 0. ? thresh : 0.5*std::sqrt(_rows+_cols+1.)*_s[0]*MathConsts::CFrealEps());
	for (j=0;j<_cols;j++) if (_s[j] > tsh) nr++;
	return nr;
}

//////////////////////////////////////////////////////////////////////////////

CFuint SVDInverter::nullity(CFreal thresh = -1.) {
	CFuint j,nn=0;
	CFreal tsh = (thresh >= 0. ? thresh : 0.5*std::sqrt(_rows+_cols+1.)*_s[0]*MathConsts::CFrealEps());
	for (j=0;j<_cols;j++) if (_s[j] <= tsh) nn++;
	return nn;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix SVDInverter::range(CFreal thresh = -1.){
	CFuint i,j,nr=0;
	CFreal tsh = (thresh >= 0. ? thresh : 0.5*std::sqrt(_rows+_cols+1.)*_s[0]*MathConsts::CFrealEps());
	RealMatrix rnge(_rows,rank(thresh));
	for (j=0;j<_cols;j++) {
		if (_s[j] > tsh) {
			for (i=0;i<_rows;i++) rnge(i,nr) = _u(i,j);
			nr++;
		}
	}
	return rnge;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix SVDInverter::nullspace(CFreal thresh = -1.){
	CFuint j,jj,nn=0;
	CFreal tsh = (thresh >= 0. ? thresh : 0.5*std::sqrt(_rows+_cols+1.)*_s[0]*MathConsts::CFrealEps());
	RealMatrix nullsp(_cols,nullity(thresh));
	for (j=0;j<_cols;j++) {
		if (_s[j] <= tsh) {
			for (jj=0;jj<_cols;jj++) nullsp(jj,nn) = _v(jj,j);
			nn++;
		}
	}
	return nullsp;
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::decompose() {
	bool flag;
	CFint i=0,its=0,j=0,jj=0,k=0,l=0,nm=0;
  CFint rows(_rows), cols(_cols); //conversion from CFuint to CFint
	CFreal anorm,c,f,g,h,s,scale,x,y,z;
	RealVector rv1(_cols);
	g = scale = anorm = 0.0;
	for (i=0;i<cols;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < rows) {
			for (k=i;k<rows;k++) scale += std::abs(_u(k,i));
			if (scale != 0.0) {
				for (k=i;k<rows;k++) {
					_u(k,i) /= scale;
					s += _u(k,i)*_u(k,i);
				}
				f=_u(i,i);
				g = -MathFunctions::changeSign(std::sqrt(s),f);
				h=f*g-s;
				_u(i,i)=f-g;
				for (j=l-1;j<cols;j++) {
					for (s=0.0,k=i;k<rows;k++) s += _u(k,i)*_u(k,j);
					f=s/h;
					for (k=i;k<rows;k++) _u(k,j) += f*_u(k,i);
				}
				for (k=i;k<rows;k++) _u(k,i) *= scale;
			}
		}
		_s[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= rows && i+1 != cols) {
			for (k=l-1;k<cols;k++) scale += std::abs(_u(i,k));
			if (scale != 0.0) {
				for (k=l-1;k<cols;k++) {
					_u(i,k) /= scale;
					s += _u(i,k)*_u(i,k);
				}
				f=_u(i,l-1);
				g = -MathFunctions::changeSign(std::sqrt(s),f);
				h=f*g-s;
				_u(i,l-1)=f-g;
				for (k=l-1;k<cols;k++) rv1[k]=_u(i,k)/h;
				for (j=l-1;j<rows;j++) {
					for (s=0.0,k=l-1;k<cols;k++) s += _u(j,k)*_u(i,k);
					for (k=l-1;k<cols;k++) _u(j,k) += s*rv1[k];
				}
				for (k=l-1;k<cols;k++) _u(i,k) *= scale;
			}
		}
		anorm=std::max(anorm,(std::abs(_s[i])+std::abs(rv1[i])));
	}
	for (i=_cols-1;i>=0;i--) {
		if (i < cols-1) {
			if (g != 0.0) {
				for (j=l;j<cols;j++) {
					_v(j,i)=(_u(i,j)/_u(i,l))/g;
				}
				for (j=l;j<cols;j++) {
					for (s=0.0,k=l;k<cols;k++) s += _u(i,k)*_v(k,j);
					for (k=l;k<cols;k++) _v(k,j) += s*_v(k,i);
				}
			}
			for (j=l;j<cols;j++) _v(i,j)=_v(j,i)=0.0;
		}
		_v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=std::min(_rows,_cols)-1;i>=0;i--) {
		l=i+1;
		g=_s[i];
		for (j=l;j<cols;j++) _u(i,j)=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<cols;j++) {
				for (s=0.0,k=l;k<rows;k++) s += _u(k,i)*_u(k,j);
				f=(s/_u(i,i))*g;
				for (k=i;k<rows;k++) _u(k,j) += f*_u(k,i);
			}
			for (j=i;j<rows;j++) _u(j,i) *= g;
		} else for (j=i;j<rows;j++) _u(j,i)=0.0;
		++_u(i,i);
	}
	for (k=_cols-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (l == 0 || std::abs(rv1[l]) <= MathConsts::CFrealEps()*anorm) {
					flag=false;
					break;
				}
				if (std::abs(_s[nm]) <= MathConsts::CFrealEps()*anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (std::abs(f) <= MathConsts::CFrealEps()*anorm) break;
					g=_s[i];
					h=pythag(f,g);
					_s[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<rows;j++) {
						y=_u(j,nm);
						z=_u(j,i);
						_u(j,nm)=y*c+z*s;
						_u(j,i)=z*c-y*s;
					}
				}
			}
			z=_s[k];
			if (l == k) {
				if (z < 0.0) {
					_s[k] = -z;
					for (j=0;j<cols;j++) _v(j,k) = -_v(j,k);
				}
				break;
			}
			if (its == 29) throw SVDException(FromHere(),"no convergence in 30 svdcmp iterations");
			x=_s[l];
			nm=k-1;
			y=_s[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+MathFunctions::changeSign(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=_s[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<cols;jj++) {
					x=_v(jj,j);
					z=_v(jj,i);
					_v(jj,j)=x*c+z*s;
					_v(jj,i)=z*c-x*s;
				}
				z=pythag(f,h);
				_s[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<rows;jj++) {
					y=_u(jj,j);
					z=_u(jj,i);
					_u(jj,j)=y*c+z*s;
					_u(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			_s[k]=x;
		}
	}
  _decomposed = true;
}

//////////////////////////////////////////////////////////////////////////////

void SVDInverter::reorder() {
	CFuint i,j,k,s,inc=1;
	CFreal sw;
	RealVector su(_rows), sv(_cols);
	do { inc *= 3; inc++; } while (inc <= _cols);
	do {
		inc /= 3;
		for (i=inc;i<_cols;i++) {
			sw = _s[i];
			for (k=0;k<_rows;k++) su[k] = _u(k,i);
			for (k=0;k<_cols;k++) sv[k] = _v(k,i);
			j = i;
			while (_s[j-inc] < sw) {
				_s[j] = _s[j-inc];
				for (k=0;k<_rows;k++) _u(k,j) = _u(k,j-inc);
				for (k=0;k<_cols;k++) _v(k,j) = _v(k,j-inc);
				j -= inc;
				if (j < inc) break;
			}
			_s[j] = sw;
			for (k=0;k<_rows;k++) _u(k,j) = su[k];
			for (k=0;k<_cols;k++) _v(k,j) = sv[k];

		}
	} while (inc > 1);
	for (k=0;k<_cols;k++) {
		s=0;
		for (i=0;i<_rows;i++) if (_u(i,k) < 0.) s++;
		for (j=0;j<_cols;j++) if (_v(j,k) < 0.) s++;
		if (s > (_rows+_cols)/2) {
			for (i=0;i<_rows;i++) _u(i,k) = -_u(i,k);
			for (j=0;j<_cols;j++) _v(j,k) = -_v(j,k);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

CFreal SVDInverter::pythag(const CFreal a, const CFreal b) {
	CFreal absa=std::abs(a), absb=std::abs(b);
	return (absa > absb ? absa*std::sqrt(1.0+ (absb/absa)*(absb/absa) ) :
		(absb == 0.0 ? 0.0 : absb*std::sqrt(1.0+ (absa/absb)*(absa/absb) )));
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
