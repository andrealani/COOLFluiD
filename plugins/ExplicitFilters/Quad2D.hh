#ifndef COOLFluiD_MathTools_Quad2D_hh
#define COOLFluiD_MathTools_Quad2D_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * Class Quad2D:
 * This class peforms simple Gaussian Quadrature integration of a
 * given function in 2 dimensions
 */
class Quad2D {  
    
public:
  
  /**
   * Integration of function func with bounds x1, x2, y1(x), y2(x)
   */
  template <class T, class T2>
  CFreal integrate(T &func, const CFreal& x1, const CFreal& x2, T2& y1, T2& y2){
    NRf1<T,T2> f1(y1,y2);
    f1.f2.func2d=&func;
    return qgaus(f1,x1,x2);
  } 
    
  /**
   * Integration of function func with constant bounds x1, x2, y1, y2
   */
  template <class T>
  CFreal integrate(T &func, const CFreal& x1, const CFreal& x2, const CFreal& y1, const CFreal& y2){
    NRf1C<T> f1(y1,y2);
    f1.f2.func2d=&func;
    return qgaus(f1,x1,x2);
  } 

private:
  template <typename T>
  struct NRf2 {
    CFreal xsav,ysav;
    // CFreal (*func2d)(const CFreal&, const CFreal&);
    T* func2d;
    CFreal operator()(const CFreal& y)
    {
        return (*func2d)(xsav,y);
    }
  };
  
  template <typename T, typename T2>
  struct NRf1 {
    // CFreal (*y1)(CFreal);
    // CFreal (*y2)(CFreal);
    T2 y1, y2;
    NRf2<T> f2;
    NRf1(T2& yy1, T2& yy2) : y1(yy1),y2(yy2) {}
    CFreal operator()(const CFreal x)
    {
        f2.xsav=x;
        return Quad2D::qgaus(f2,y1(x),y2(x));
    }
  };
  
  template <typename T>
  struct NRf1C {
    NRf2<T> f2;
    CFreal y1, y2;
    NRf1C(CFreal yy1, CFreal yy2) : y1(yy1),y2(yy2) {}
    CFreal operator()(const CFreal x)
    {
        f2.xsav=x;
        return Quad2D::qgaus(f2,y1,y2);
    }
  };
    
  static const CFreal x[];
  static const CFreal w[];
  
  template <class T>
  static CFreal qgaus(T &func, const CFreal a, const CFreal b)
  {
    CFreal xm=0.5*(b+a);
    CFreal xr=0.5*(b-a);
    CFreal s=0;
    for (CFuint j=0;j<5;j++) {
        CFreal dx=xr*x[j];
        s += w[j]*(func(xm+dx)+func(xm-dx));
    }
    return s *= xr;
  }
};

const CFreal Quad2D::x[]={0.1488743389816312,0.4333953941292472,
0.6794095682990244,0.8650633666889845,0.9739065285171717};
const CFreal Quad2D::w[]={0.2955242247147529,0.2692667193099963,
0.2190863625159821,0.1494513491505806,0.0666713443086881};

//////////////////////////////////////////////////////////////////////////////

  } // end namespace MathTools

} // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_Quad2D_hh
