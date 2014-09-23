// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_FindMinimum_hh
#define COOLFluiD_MathTools_FindMinimum_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * BracketMethod:
 * This class is the basis class to find the minimum of a given function.
 * This basis class allows to bracket the minimum and assert that the 
 * minimum exists.
 * 
 * @author Willem Deconinck
 * Based on Numerical Recipes v3.02
 * 
 * Use:
 * @verbatim
 * class Function
 * {
 *   operator() (const CFreal& x) { return x*x;}
 * };
 * Function func;
 * BracketMethod bracket;
 * bracket.bracket(a,b,func);
 * @endverbatim
 */
class BracketMethod {
public: 
  
  /**
   * Default constructor
   */
  BracketMethod() : m_maxStep(MathConsts::CFrealMax()) {}
  
  /**
   * Bracket the minimum in a given function, 
   * starting with 2 values pointing downwards 
   * to the minimum
   * @param a A first value 
   * @param b A second value
   * @param func The function to bracket
   */
	template <class T>
	void bracket(const CFreal a, const CFreal b, T& func)
	{
		const CFreal GLIMIT=100.0,TINY=1.0e-20;
    const CFreal GOLD=1.618034;
		ax=a; bx=b;
		CFreal fu;
		fa=func(ax);
		fb=func(bx);
		if (fb > fa) {
			swap(ax,bx);
			swap(fb,fa);
		}
		cx=bx+std::min(GOLD*(bx-ax),m_maxStep);
		fc=func(cx);
		while (fb > fc) {
			CFreal r=(bx-ax)*(fb-fc);
			CFreal q=(bx-cx)*(fb-fa);
			CFreal u=bx-((bx-cx)*q-(bx-ax)*r)/
				(2.0*MathTools::MathFunctions::changeSign(std::max(std::abs(q-r),TINY),q-r));
			CFreal ulim=bx+GLIMIT*(cx-bx);
			if ((bx-u)*(u-cx) > 0.0) {
				fu=func(u);
				if (fu < fc) {
					ax=bx;
					bx=u;
					fa=fb;
					fb=fu;
					return;
				} else if (fu > fb) {
					cx=u;
					fc=fu;
					return;
				}
				u=cx+std::min(GOLD*(cx-bx),m_maxStep);
				fu=func(u);
			} else if ((cx-u)*(u-ulim) > 0.0) {
				fu=func(u);
				if (fu < fc) {
					shft3(bx,cx,u,u+std::min(GOLD*(u-cx),m_maxStep));
					shft3(fb,fc,fu,func(u));
				}
			} else if ((u-ulim)*(ulim-cx) >= 0.0) {
				u=ulim;
				fu=func(u);
			} else {
				u=cx+std::min(GOLD*(cx-bx),m_maxStep);
				fu=func(u);
			}
			shft3(ax,bx,cx,u);
			shft3(fa,fb,fc,fu);
		}
	}
	
	/**
   * Sets the maximum step taken to bracket the minimum.
   * This is purely optional
   * @param step the maximum step
   */
  inline void setMaxStep(const CFreal& step)
  {
    m_maxStep = step;
  }
  
protected:
  CFreal m_maxStep;
 	CFreal ax,bx,cx,fa,fb,fc;
 	  
	inline void shft2(CFreal &a, CFreal &b, const CFreal c)
	{
		a=b;
		b=c;
	}
	inline void shft3(CFreal &a, CFreal &b, CFreal &c, const CFreal d)
	{
		a=b;
		b=c;
		c=d;
	}
	inline void mov3(CFreal &a, CFreal &b, CFreal &c, const CFreal d, const CFreal e, const CFreal f)
	{
		a=d; b=e; c=f;
	}
	
  inline void swap(CFreal &a, CFreal &b)
  {
    CFreal dum=a; a=b; b=dum;
  }
	
};

//////////////////////////////////////////////////////////////////////////////

/**
 * Golden:
 * This class is finds the minimum of a given function, 
 * using the golden method.
 * 
 * @author Willem Deconinck
 * Based on Numerical Recipes v3.02
 * 
 * Use:
 * @verbatim
 * class Function
 * {
 *   CFreal operator() (const CFreal& x) { return x*x;}
 * };
 * Function func;
 * Golden golden;
 * golden.bracket(a,b,func);
 * CFreal xmin = golden.minimize(func);
 * @endverbatim
 */
class Golden : public BracketMethod {
private:
	CFreal xmin,fmin;
	const CFreal tol;
public:
  
  /**
   * Default constructor
   */
	Golden(const CFreal toll=3.0e-8) : BracketMethod(), tol(toll) {}
	
	/**
   * Minimizes the given function
   * @param func The given function
   * @return The minimum of the function
   */
	template <class T>
	CFreal minimize(T &func)
	{
		const CFreal R=0.61803399,C=1.0-R;
		CFreal x1,x2;
		CFreal x0=ax;
		CFreal x3=cx;
		if (std::abs(cx-bx) > std::abs(bx-ax)) {
			x1=bx;
			x2=bx+C*(cx-bx);
		} else {
			x2=bx;
			x1=bx-C*(bx-ax);
		}
		CFreal f1=func(x1);
		CFreal f2=func(x2);
		while (std::abs(x3-x0) > tol*(std::abs(x1)+std::abs(x2))) {
			if (f2 < f1) {
				shft3(x0,x1,x2,R*x2+C*x3);
				shft2(f1,f2,func(x2));
			} else {
				shft3(x3,x2,x1,R*x1+C*x0);
				shft2(f2,f1,func(x1));
			}
		}
		if (f1 < f2) {
			xmin=x1;
			fmin=f1;
		} else {
			xmin=x2;
			fmin=f2;
		}
		return xmin;
	}
};

//////////////////////////////////////////////////////////////////////////////

/**
 * Brent:
 * This class is finds the minimum of a given function, 
 * using the Brent method.
 * 
 * @author Willem Deconinck
 * Based on Numerical Recipes v3.02
 * 
 * Use:
 * @verbatim
 * class Function
 * {
 *   CFreal operator() (const CFreal& x) { return x*x;}
 * };
 * Function func;
 * Brent brent;
 * brent.bracket(a,b,func);
 * CFreal xmin = brent.minimize(func);
 * @endverbatim
 */
class Brent : public BracketMethod {
private:
	CFreal xmin,fmin;
	const CFreal tol;
public:
  
  /**
   * Default constructor
   */
  Brent(const CFreal toll=3.0e-8) : BracketMethod(), tol(toll) {}
  
  /**
   * Minimizes the given function
   * @param func The given function
   * @return The minimum of the function
   */
	template <class T>
	CFreal minimize(T &func)
	{
		const CFuint ITMAX=100;
		const CFreal CGOLD=0.3819660;
		const CFreal ZEPS=MathTools::MathConsts::CFrealEps()*1.0e-3;
		CFreal a,b,d=0.0,etemp,fu,fv,fw,fx;
		CFreal p,q,r,tol1,tol2,u,v,w,x,xm;
		CFreal e=0.0;
	
		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=func(x);
		for (CFuint iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*std::abs(x)+ZEPS);
			if (std::abs(x-xm) <= (tol2-0.5*(b-a))) {
				fmin=fx;
				return xmin=x;
			}
			if (std::abs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=std::abs(q);
				etemp=e;
				e=d;
				if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q*(a-x)
						|| p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=MathTools::MathFunctions::changeSign(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(std::abs(d) >= tol1 ? x+d : x+MathTools::MathFunctions::changeSign(tol1,d));
			fu=func(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				shft3(v,w,x,u);
				shft3(fv,fw,fx,fu);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		throw("Too many iterations in brent");
	}
};

//////////////////////////////////////////////////////////////////////////////

/**
 * Dbrent:
 * This class is finds the minimum of a given function, 
 * using the Brent method, and making use of the first
 * derivative of the function.
 * The given function-object requires a subroutine df(x)
 * 
 * @author Willem Deconinck
 * Based on Numerical Recipes v3.02
 * 
 * Use:
 * @verbatim
 * class Function
 * {
 *   CFreal operator() (const CFreal& x) { return x*x;}
 *   CFreal df(const CFreal& x) { return 2*x;}
 * };
 * Function func;
 * Dbrent dbrent;
 * dbrent.bracket(a,b,func);
 * CFreal xmin = dbrent.minimize(func);
 * @endverbatim
 */
class Dbrent : public BracketMethod {
private:
	CFreal xmin,fmin;
	const CFreal tol;
public:
  
  /**
   * Default constructor
   */
  Dbrent(const CFreal toll=3.0e-8) : BracketMethod(), tol(toll) {}
  
  /**
   * Minimizes the given function
   * @param funcd The given function
   * @return The minimum of the function
   */
	template <class T>
	CFreal minimize(T &funcd)
	{
		const CFuint ITMAX=100;
		const CFreal ZEPS=MathTools::MathConsts::CFrealEps()*1.0e-3;
		bool ok1,ok2;
		CFreal a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
		CFreal fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	
		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=funcd(x);
		dw=dv=dx=funcd.df(x);
		for (CFuint iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol1=tol*std::abs(x)+ZEPS;
			tol2=2.0*tol1;
			if (std::abs(x-xm) <= (tol2-0.5*(b-a))) {
				fmin=fx;
				return xmin=x;
			}
			if (std::abs(e) > tol1) {
				d1=2.0*(b-a);
				d2=d1;
				if (dw != dx) d1=(w-x)*dx/(dx-dw);
				if (dv != dx) d2=(v-x)*dx/(dx-dv);
				u1=x+d1;
				u2=x+d2;
				ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
				ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				olde=e;
				e=d;
				if (ok1 || ok2) {
					if (ok1 && ok2)
						d=(std::abs(d1) < std::abs(d2) ? d1 : d2);
					else if (ok1)
						d=d1;
					else
						d=d2;
					if (std::abs(d) <= std::abs(0.5*olde)) {
						u=x+d;
						if (u-a < tol2 || b-u < tol2)
							d=MathTools::MathFunctions::changeSign(tol1,xm-x);
					} else {
						d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
					}
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
			if (std::abs(d) >= tol1) {
				u=x+d;
				fu=funcd(u);
			} else {
				u=x+MathTools::MathFunctions::changeSign(tol1,d);
				fu=funcd(u);
				if (fu > fx) {
					fmin=fx;
					return xmin=x;
				}
			}
			du=funcd.df(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				mov3(v,fv,dv,w,fw,dw);
				mov3(w,fw,dw,x,fx,dx);
				mov3(x,fx,dx,u,fu,du);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					mov3(v,fv,dv,w,fw,dw);
					mov3(w,fw,dw,u,fu,du);
				} else if (fu < fv || v == x || v == w) {
					mov3(v,fv,dv,u,fu,du);
				}
			}
		}
		throw("Too many iterations in routine dbrent");
	}
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_FindMinimum_hh

