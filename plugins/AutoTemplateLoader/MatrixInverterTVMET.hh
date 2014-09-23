#ifndef COOLFluiD_MathTools_MatrixInverterTVMET_hh
#define COOLFluiD_MathTools_MatrixInverterTVMET_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 

  namespace MathTools { 

//////////////////////////////////////////////////////////////////////////////

/// This class inverts a generic matrix.
/// It takes a parameter with the size of the square matrix.
/// It can be specialised for small sizes.
///
/// @author Andrea Lani
/// @author Tiago Quintino
template < unsigned int SIZE >
struct MatrixInverterTVMET {

  /// Invert the given matrix a and put the result in x
  void invert (const tvmet::Matrix<CFreal,SIZE,SIZE>& a, tvmet::Matrix<CFreal,SIZE,SIZE>& x)
  {
    throw std::string ("MatrixInverterTVMET not implemented\n");
  }

private: // data

}; // class MatrixInverterTVMET

//////////////////////////////////////////////////////////////////////////////

/// Specialization of matrix inversion for SIZE 2
/// @author Tiago Quintino
template <> struct MatrixInverterTVMET <2> {

  /// Invert the given matrix a and put the result in x
  void invert (const tvmet::Matrix<CFreal,2,2>& a, tvmet::Matrix<CFreal,2,2>& x)
  {
    cf_assert(a.rows() == 2);
    cf_assert(a.cols() == 2);
    cf_assert(x.rows() == 2);
    cf_assert(x.cols() == 2);

    const CFreal det = a(0,0)*a(1,0) - a(0,1)*a(1,1); cf_assert(MathChecks::isNotZero(det));
    const CFreal invDet = 1./det;

    x(0,0) =  a(1,1)*invDet;
    x(0,1) = -a(0,1)*invDet;
    x(1,0) = -a(1,0)*invDet;
    x(1,1) =  a(0,0)*invDet;
  }

}; // class MatrixInverterTVMET<2>

//////////////////////////////////////////////////////////////////////////////

/// Specialization of matrix inversion for SIZE 3
/// @author Tiago Quintino
template <> struct MatrixInverterTVMET <3> {

  /// Invert the given matrix a and put the result in x
  void invert (const tvmet::Matrix<CFreal,3,3>& a, tvmet::Matrix<CFreal,3,3>& x)
  {
    cf_assert(a.rows() == 3);
    cf_assert(a.cols() == 3);
    cf_assert(x.rows() == 3);
    cf_assert(x.cols() == 3);

    const CFreal det = a(0,0)*(a(1,1)*a(2,2) - a(1,2)*a(2,1)) -
                       a(0,1)*(a(1,0)*a(2,2) - a(1,2)*a(2,0)) +
                       a(0,2)*(a(1,0)*a(2,1) - a(1,1)*a(2,0));  cf_assert(MathChecks::isNotZero(det));

    const CFreal invDet = 1./det;

    x(0,0) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1))*invDet;
    x(0,1) = -(a(0,1)*a(2,2) - a(0,2)*a(2,1))*invDet;
    x(0,2) =  (a(0,1)*a(1,2) - a(1,1)*a(0,2))*invDet;
    x(1,0) = -(a(1,0)*a(2,2) - a(1,2)*a(2,0))*invDet;
    x(1,1) =  (a(0,0)*a(2,2) - a(0,2)*a(2,0))*invDet;
    x(1,2) = -(a(0,0)*a(1,2) - a(0,2)*a(1,0))*invDet;
    x(2,0) =  (a(1,0)*a(2,1) - a(1,1)*a(2,0))*invDet;
    x(2,1) = -(a(0,0)*a(2,1) - a(0,1)*a(2,0))*invDet;
    x(2,2) =  (a(0,0)*a(1,1) - a(0,1)*a(1,0))*invDet;
  }

}; // class MatrixInverterTVMET<3>

//////////////////////////////////////////////////////////////////////////////

/// Specialization of matrix inversion for SIZE 4
/// @author Tiago Quintino
template <> struct MatrixInverterTVMET <4> {

  /// Invert the given matrix a and put the result in x
  void invert (const tvmet::Matrix<CFreal,4,4>& a, tvmet::Matrix<CFreal,4,4>& x)
  {
    cf_assert(a.rows() == 4);
    cf_assert(a.cols() == 4);
    cf_assert(x.rows() == 4);
    cf_assert(x.cols() == 4);

    const CFreal t14 = a(0,0)*a(1,1) ;
    const CFreal t15 = a(2,2)*a(3,3) ;
    const CFreal t17 = a(2,3)*a(3,2) ;
    const CFreal t19 = a(0,0)*a(2,1) ;
    const CFreal t20 = a(1,2)*a(3,3) ;
    const CFreal t22 = a(1,3)*a(3,2) ;
    const CFreal t24 = a(0,0)*a(3,1) ;
    const CFreal t25 = a(1,2)*a(2,3) ;
    const CFreal t27 = a(1,3)*a(2,2) ;
    const CFreal t29 = a(1,0)*a(0,1) ;
    const CFreal t32 = a(1,0)*a(2,1) ;
    const CFreal t33 = a(0,2)*a(3,3) ;
    const CFreal t35 = a(0,3)*a(3,2) ;
    const CFreal t37 = a(1,0)*a(3,1) ;
    const CFreal t38 = a(0,2)*a(2,3) ;
    const CFreal t40 = a(0,3)*a(2,2) ;
    const CFreal t42 = t14*t15-t14*t17-t19*t20+t19*t22+t24*t25-t24*t27-t29*t15+t29*t17+t32*t33-t32*t35-t37*t38+t37*t40;
    const CFreal t43 = a(2,0)*a(0,1) ;
    const CFreal t46 = a(2,0)*a(1,1) ;
    const CFreal t49 = a(2,0)*a(3,1) ;
    const CFreal t50 = a(0,2)*a(1,3) ;
    const CFreal t52 = a(0,3)*a(1,2) ;
    const CFreal t54 = a(3,0)*a(0,1) ;
    const CFreal t57 = a(3,0)*a(1,1) ;
    const CFreal t60 = a(3,0)*a(2,1) ;
    const CFreal t63 = t43*t20-t43*t22-t46*t33+t46*t35+t49*t50-t49*t52-t54*t25+t54*t27+t57*t38-t57*t40-t60*t50+t60*t52;

    const CFreal ddet = t42+t63 ;
    cf_assert(MathChecks::isNotZero(ddet));
    const CFreal deter = 1. / ddet ;

    const CFreal t71  = a(2,1)*a(0,2) ;
    const CFreal t73  = a(2,1)*a(0,3) ;
    const CFreal t75  = a(3,1)*a(0,2) ;
    const CFreal t77  = a(3,1)*a(0,3) ;
    const CFreal t81  = a(0,1)*a(1,2) ;
    const CFreal t83  = a(0,1)*a(1,3) ;
    const CFreal t85  = a(1,1)*a(0,2) ;
    const CFreal t87  = a(1,1)*a(0,3) ;
    const CFreal t119 = a(2,0)*a(0,2) ;
    const CFreal t121 = a(2,0)*a(0,3) ;
    const CFreal t123 = a(3,0)*a(0,2) ;
    const CFreal t125 = a(3,0)*a(0,3) ;
    const CFreal t129 = a(0,0)*a(1,2) ;
    const CFreal t131 = a(0,0)*a(1,3) ;
    const CFreal t133 = a(1,0)*a(0,2) ;
    const CFreal t135 = a(1,0)*a(0,3) ;

    x(0,0) = (a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)*a(3,3)
          +  a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)*a(1,3)*a(2,2))*deter ;
    x(0,1) = -(a(0,1)*a(2,2)*a(3,3)-a(0,1)*a(2,3)*a(3,2)-t71*a(3,3)+t73*a(3,2)
        +   t75*a(2,3)-t77*a(2,2))*deter ;
    x(0,2) = (t81*a(3,3)-t83*a(3,2)-t85*a(3,3)+t87*a(3,2)+t75*a(1,3)-t77*a(1,2))*deter ;
    x(0,3) = -(t81*a(2,3)-t83*a(2,2)-t85*a(2,3)+t87*a(2,2)+t71*a(1,3)-t73*a(1,2))*deter ;

    x(1,0) = (-a(1,0)*a(2,2)*a(3,3)+a(1,0)*a(2,3)*a(3,2)+a(2,0)*a(1,2)*a(3,3)
                -   a(2,0)*a(1,3)*a(3,2)-a(3,0)*a(1,2)*a(2,3)+a(3,0)*a(1,3)*a(2,2))*deter ;
    x(1,1) = (a(0,0)*a(2,2)*a(3,3)-a(0,0)*a(2,3)*a(3,2)-t119*a(3,3)+t121*a(3,2)
                +  t123*a(2,3)-t125*a(2,2))*deter ;
    x(1,2) = -(t129*a(3,3)-t131*a(3,2)-t133*a(3,3)+t135*a(3,2)+t123*a(1,3)-t125*a(1,2))*deter ;
    x(1,3) = (t129*a(2,3)-t131*a(2,2)-t133*a(2,3)+t135*a(2,2)+t119*a(1,3)-t121*a(1,2))*deter ;

    x(2,0) = -(-t32*a(3,3)+t37*a(2,3)+t46*a(3,3)-t49*a(1,3)-t57*a(2,3)+t60*a(1,3))*deter ;
    x(2,1) = -(t19*a(3,3)-t24*a(2,3)-t43*a(3,3)+t121*a(3,1)+t54*a(2,3)-t125*a(2,1))*deter ;
    x(2,2) = (t14*a(3,3)-t24*a(1,3)-t29*a(3,3)+t135*a(3,1)+t54*a(1,3)-t125*a(1,1))*deter ;
    x(2,3) = -(t14*a(2,3)-t19*a(1,3)-t29*a(2,3)+t135*a(2,1)+t43*a(1,3)-t121*a(1,1))*deter ;

    x(3,0) = -(t32*a(3,2)-t37*a(2,2)-t46*a(3,2)+t49*a(1,2)+t57*a(2,2)-t60*a(1,2))*deter ;
    x(3,1) = (t19*a(3,2)-t24*a(2,2)-t43*a(3,2)+t119*a(3,1)+t54*a(2,2)-t123*a(2,1))*deter ;
    x(3,2) = -(t14*a(3,2)-t24*a(1,2)-t29*a(3,2)+t133*a(3,1)+t54*a(1,2)-t123*a(1,1))*deter ;
    x(3,3) = (t14*a(2,2)-t19*a(1,2)-t29*a(2,2)+t133*a(2,1)+t43*a(1,2)-t119*a(1,1))*deter ;
  }

}; // class MatrixInverterTVMET<4>

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD 

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_MatrixInverterTVMET_hh
