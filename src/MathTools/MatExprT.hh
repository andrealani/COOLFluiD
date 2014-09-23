// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MatExprT_hh
#define COOLFluiD_MathTools_MatExprT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/ExprT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
  
//////////////////////////////////////////////////////////////////////////////

/**
 * Definition of a wrapper base expression class MatExprT which deruves from 
 * ExprT. The expression template accepts two template parameters:
 * 1- the derived expression (array)
 * 2- a tuple argument storing all arguments of the derived expression: those 
 *    arguments cannot be provided directly by DATA (via typedefs or enum) because
 *    at the moment of instantiation, DATA is still an incomplete type, being 
 *    deriving from ExprT itself.
 *
 * @author Andrea Lani
 */
template <typename DATA, typename ARG>
class MatExprT : public ExprT<DATA, ARG> {
public:
  enum {SIZE2=ARG::SIZE2};
  
  /// Constructor
  HHOST_DEV MatExprT(DATA* data) : ExprT<DATA, ARG>(data) {}
  
  /// Destructor
  HHOST_DEV ~MatExprT(){} 
  
  /// Number of columns
  HHOST_DEV size_t nbCols() const {return this->getData()->nbCols();}
  
  /// Number of rows
  HHOST_DEV size_t nbRows() const {return this->getData()->nbRows();}

  /// Accessor to individual entry
  HHOST_DEV typename ARG::TYPE at(size_t i, size_t j) const {return this->getData()->at(i,j);} 
  
};
  
//////////////////////////////////////////////////////////////////////////////

/**
 * Definition of an expression template class for basic binary operations.
 * A macro is used to save tedious code duplication.
 * The expression template accepts two parameters:
 * 1. first operand expression type
 * 2. second operand expression type.
 * The constructor accepts ET objects (corresponding to the two operands of 
 * the expression) as arguments.
 *
 * @author Andrea Lani
 */
#define EMET_BINARY(__OpName__,__op__) \
template <typename V1, typename V2>			\
class __OpName__ : public MatExprT< __OpName__<V1,V2>, ETPLMAT(V1,V2)> { \
public:						\
  HHOST_DEV __OpName__ (MEETYPE(V1) v1, MEETYPE(V2) v2) :	    \
    MatExprT<__OpName__<V1,V2>, ETPLMAT(V1,V2)>(this), e1(v1), e2(v2) {}	\
    									\
    HHOST_DEV ETPLTYPE(V1) at(size_t i, size_t j) const {return e1.at(i,j) __op__ e2.at(i,j);} \
    									\
    HHOST_DEV size_t size()   const {return e1.size();}		\
    HHOST_DEV size_t nbCols() const {return e1.nbCols();}        \
    HHOST_DEV size_t nbRows() const {return e1.nbRows();}        \
									\
    MEETYPE(V1) e1;							\
    MEETYPE(V2) e2;							\
};
  
EMET_BINARY(MAddT,+)
EMET_BINARY(MSubT,-)
EMET_BINARY(MDivT,/)

#undef EMET_BINARY
  
//////////////////////////////////////////////////////////////////////////////

/**
 * Definition of an expression template class for basic binary operations.
 * A macro is used to save tedious code duplication.
 * The expression template accepts one parameter, corresponding to the 
 * first operand expression type.
 * The constructor accepts one ET object and a constant value (corresponding 
 * to the two operands of the expression) as arguments.
 *
 * @author Andrea Lani
 */
#define EMET_BINARY_C(__a1__,__a2__,__e1__,__e2__,__e__,__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public MatExprT< __OpName__<V>, ETPL(V)> { \
public:						\
  HHOST_DEV __OpName__ (__a1__ v1, __a2__ v2) :	    \
    MatExprT<__OpName__<V>, ETPL(V)>(this), e1(v1), e2(v2) {}	\
    									\
    HHOST_DEV ETPLTYPE(V) at(size_t i, size_t j) const {return __e1__ __op__ __e2__;} \
    									\
    HHOST_DEV size_t size()   const {return __e__ .size();}		\
    HHOST_DEV size_t nbCols() const {return __e__ .nbCols();}		\
    HHOST_DEV size_t nbRows() const {return __e__ .nbRows();}		\
									\
    __a1__ e1;								\
    __a2__ e2;								\
};
  
EMET_BINARY_C(MEETYPE(V), ETPLTYPE(V), e1.at(i,j), e2, e1, MAddTv1,+)
EMET_BINARY_C(MEETYPE(V), ETPLTYPE(V), e1.at(i,j), e2, e1, MSubTv1,-)
EMET_BINARY_C(MEETYPE(V), ETPLTYPE(V), e1.at(i,j), e2, e1, MMulTv1,*)
EMET_BINARY_C(MEETYPE(V), ETPLTYPE(V), e1.at(i,j), e2, e1, MDivTv1,/)

EMET_BINARY_C(ETPLTYPE(V), MEETYPE(V), e1, e2.at(i,j), e2, MAddTv2,+)
EMET_BINARY_C(ETPLTYPE(V), MEETYPE(V), e1, e2.at(i,j), e2, MSubTv2,-)
EMET_BINARY_C(ETPLTYPE(V), MEETYPE(V), e1, e2.at(i,j), e2, MMulTv2,*)
EMET_BINARY_C(ETPLTYPE(V), MEETYPE(V), e1, e2.at(i,j), e2, MDivTv2,/)

#undef EMET_BINARY_C
  
//////////////////////////////////////////////////////////////////////////////

/**
 * Inline function (the keyword must be put otherwise the compiler will not inline it,
 * as shown by the profiler) providing the overloading of the basic binary operators
 * between two expression template objects.
 *
 * @param v1 expression template object (first operand)
 * @param v2 expression template object (second operand)
 */

//////////////////////////////////////////////////////////////////////////////

#define EMET_BINARY_OP(__t1__,__t2__,__a1__,__a2__,__OpName__,__op__)	\
template __t1__							\
HHOST_DEV  inline __OpName__ __t2__ operator __op__ (__a1__ v1, __a2__ v2) \
{									\
  return __OpName__ __t2__(v1,v2);						\
}

/// binary operator overloading taking two expressions
EMET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, MEETYPE(V1), MEETYPE(V2), MAddT, +)
EMET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, MEETYPE(V1), MEETYPE(V2), MSubT, -)
EMET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, MEETYPE(V1), MEETYPE(V2), MDivT, /)

/// binary operator overloading taking one expression and one constant value
EMET_BINARY_OP(<typename V>, <V>, MEETYPE(V), ETPLTYPE(V), MAddTv1, +)
EMET_BINARY_OP(<typename V>, <V>, MEETYPE(V), ETPLTYPE(V), MSubTv1, -)
EMET_BINARY_OP(<typename V>, <V>, MEETYPE(V), ETPLTYPE(V), MMulTv1, *)
EMET_BINARY_OP(<typename V>, <V>, MEETYPE(V), ETPLTYPE(V), MDivTv1, /)
  
/// binary operator overloading taking one constant value and one expression
EMET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), MEETYPE(V), MAddTv2, +)
EMET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), MEETYPE(V), MSubTv2, -)
EMET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), MEETYPE(V), MMulTv2, *)
EMET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), MEETYPE(V), MDivTv2, /)
  
#undef EMET_BINARY_OP 

//////////////////////////////////////////////////////////////////////////////

/**
 * Definition of an expression template class for basic binary operations.
 * A macro is used to save tedious code duplication.
 * The expression template accepts one parameter, corresponding to the 
 * first operand expression type.
 * The constructor accepts one ET object and a constant value (corresponding 
 * to the two operands of the expression) as arguments.
 *
 * @author Andrea Lani
 */
#define EMET_UNARY(__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public MatExprT< __OpName__<V>, ETPL(V)> { \
public:						\
  HHOST_DEV __OpName__ (MEETYPE(V) v1) :	    \
    MatExprT<__OpName__<V>, ETPL(V)>(this), e1(v1) {}	\
    									\
    HHOST_DEV ETPLTYPE(V) at(size_t i, size_t j) const {return __op__(e1.at(i,j));} \
    									\
    HHOST_DEV size_t size() const {return e1.size();}			\
    HHOST_DEV size_t nbCols() const {return e1.nbCols();}        \
    HHOST_DEV size_t nbRows() const {return e1.nbRows();}		\
    									\
    MEETYPE(V) e1;							\
};

  EMET_UNARY(MCosT,     std::cos)
  EMET_UNARY(MSinT,     std::sin)
  EMET_UNARY(MTanT,     std::tan)
  EMET_UNARY(MAcosT,    std::acos)
  EMET_UNARY(MAsinT,    std::asin)
  EMET_UNARY(MAtanT,    std::atan)
  EMET_UNARY(MCoshT,    std::cosh)
  EMET_UNARY(MSinhT,    std::sinh)
  EMET_UNARY(MTanhT,    std::tanh)
  EMET_UNARY(MLogT,     std::log)
  EMET_UNARY(MLog10T,   std::log10)
  EMET_UNARY(MExpT,     std::exp)
  EMET_UNARY(MSqrtT,    std::sqrt)
  EMET_UNARY(MAbsT,     std::abs)
  EMET_UNARY(MMinusOpT, -)

#undef EMET_UNARY
  
//////////////////////////////////////////////////////////////////////////////
  
#define EMET_UNARY_OP(__OpName__,__op__)		\
template <typename V>						\
HHOST_DEV inline __OpName__<V> __op__ (MEETYPE(V) v1)		\
{									\
  return __OpName__<V>(v1);						\
}
    
  /// unary operator overloading taking one expression
  /// (handle will care: they occasionally might not work on GPU)
  EMET_UNARY_OP(MSinT,     sin)
  EMET_UNARY_OP(MCosT,     cos)
  EMET_UNARY_OP(MTanT,     tan) 
  EMET_UNARY_OP(MAsinT,    asin)
  EMET_UNARY_OP(MAcosT,    acos)
  EMET_UNARY_OP(MAtanT,    atan) 
  EMET_UNARY_OP(MSinhT,    sinh)
  EMET_UNARY_OP(MCoshT,    cosh)
  EMET_UNARY_OP(MTanhT,    tanh) 
  EMET_UNARY_OP(MLogT,     log)
  EMET_UNARY_OP(MLog10T,   log10) 
  EMET_UNARY_OP(MExpT,     exp)
  EMET_UNARY_OP(MSqrtT,    sqrt) 
  EMET_UNARY_OP(MAbsT,     abs)
  EMET_UNARY_OP(MMinusOpT, operator-)

#undef EMET_UNARY_OP 

//////////////////////////////////////////////////////////////////////////////
 
/**
 * Expression class for multiplication between two CFMat-based expressions 
 * (V1 nd V2). The constructor accepts expressions V1 and V2 as arguments.
 *
 * @author Andrea Lani
 */  
template <typename V1, typename V2>
class MMulT :public MatExprT< MMulT<V1,V2>, ETPLMATMULT(V1,V2)> {
public:	
  
  /// Constructor
  HHOST_DEV MMulT (MEETYPE(V1) v1, MEETYPE(V2) v2) :
    MatExprT< MMulT<V1,V2>, ETPLMATMULT(V1,V2)>(this), ex1(v1), ex2(v2) {}
    
  /// Accessor to individual entry
  HHOST_DEV ETPLTYPE(V2) at(size_t i, size_t j) const  
  { 
    ETPLTYPE(V2) res = ETPLTYPE(V2)();
    const size_t nbCols = ex1.nbCols();
    for (size_t k = 0; k < nbCols; ++k) {
      res += ex1.at(i,k)*ex2.at(k,j);
    }
    return res;
  }
  
  /// Size of the array
  HHOST_DEV size_t size() const {return ex1.nbRows()*ex2.nbCols();}
  
  /// Number of columns
  HHOST_DEV size_t nbCols() const {return ex2.nbCols();}
  
  /// Number of rows
  HHOST_DEV size_t nbRows() const {return ex1.nbRows();}
  
private:
  
  MEETYPE(V1) ex1;			
  MEETYPE(V2) ex2;			
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename V1, typename V2>
HHOST_DEV inline MMulT<V1,V2> operator* (MEETYPE(V1) v1, MEETYPE(V2) v2)
{
  return MMulT<V1,V2>(v1,v2);
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Expression class for multiplication between a vector-based (V1) and a CFMat-based (V2)
 * expressions. It can only be used with squared matrices. V1 is intended as a diagonal 
 * matrix. The constructor accepts expressions V1 and V2 as arguments.
 *
 * @author Andrea Lani
 */  
template <typename V1, typename V2>
class VMMulT :public MatExprT< VMMulT<V1,V2>, ETPLMATMULT(V1,V2)> {
public:	
  
  /// Constructor
  HHOST_DEV VMMulT (EETYPE(V1) v1, MEETYPE(V2) v2) :
    MatExprT< VMMulT<V1,V2>, ETPLMATMULT(V1,V2)>(this), ex1(v1), ex2(v2) {}
  
  /// Accessor to individual entry
  HHOST_DEV ETPLTYPE(V2) at(size_t i, size_t j) const {return ex1.at(i)*ex2.at(i,j);}
  
  /// Size of the array
  HHOST_DEV size_t size() const {return ex2.size();}
  
  /// Number of columns
  HHOST_DEV size_t nbCols() const {return ex2.nbCols();}
  
  /// Number of rows
  HHOST_DEV size_t nbRows() const {return ex2.nbRows();}
 
private:
  
  EETYPE(V1)   ex1;			
  MEETYPE(V2)  ex2;			
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename V1, typename V2>
HHOST_DEV inline VMMulT<V1,V2> operator* (EETYPE(V1) v1, MEETYPE(V2) v2)
{
  return VMMulT<V1,V2>(v1,v2);
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Expression class for multiplication between a Mat-based (V1) and a vector-based (V2)
 * expressions. The constructor accepts expressions V1 and V2 as arguments.
 *
 * @author Andrea Lani
 */ 
template <typename V1, typename V2>
class MVMulVT :public ExprT< MVMulVT<V1,V2>, ETPLVEC(V1,V2)> {
public:	
  
  /// Constructor
  HHOST_DEV MVMulVT (MEETYPE(V1) v1, EETYPE(V2) v2) :
    ExprT< MVMulVT<V1,V2>, ETPLVEC(V1,V2)>(this), ex1(v1), ex2(v2) {}
    
  /// Accessor to individual entry
  HHOST_DEV ETPLTYPE(V2) at(size_t i) const  
  { 
    ETPLTYPE(V2) res = ETPLTYPE(V2)();
    for (size_t k = 0; k < ex2.size(); ++k) {
      res += ex1.at(i,k) * ex2.at(k);
    }
    return res;
  }
  
  /// Size of the array
  HHOST_DEV size_t size() const {return ex2.size();}

private:
  
  MEETYPE(V1) ex1;			
  EETYPE(V2)  ex2;			
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename V1, typename V2>
HHOST_DEV inline MVMulVT<V1,V2> operator* (MEETYPE(V1) v1, EETYPE(V2) v2)
{
  return MVMulVT<V1,V2>(v1,v2);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
