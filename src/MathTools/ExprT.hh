// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_ExprT_hh
#define COOLFluiD_MathTools_ExprT_hh

//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "MathTools/MacrosET.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
  
//////////////////////////////////////////////////////////////////////////////

template <typename T, int N, int M>
struct ETuple {
  typedef T TYPE;
  enum {SIZE1=N};
  enum {SIZE2=M};
};
  
//////////////////////////////////////////////////////////////////////////////

/**
 * Definition of a wrapper base expression class ExprT, the closure objects
 * (derived from ExprT using the CRTP technique), some partial or explicit 
 * specializations of the closure objects to treat special cases and 
 * definition of the overloading functions accepting expressions or constants
 * as arguments. The expression template accepts two template parameters:
 * 1- the derived expression (array)
 * 2- a tuple argument storing all arguments of the derived expression: those 
 *    arguments cannot be provided directly by DATA (via typedefs or enum) because
 *    at the moment of instantiation, DATA is still an incomplete type, being 
 *    deriving from ExprT itself.
 *
 * @author Andrea Lani
 */
template <typename DATA, typename ARG>
class ExprT {
public:
  typedef DATA PTR;
  typedef ARG TUPLE;
  enum {SIZE1=ARG::SIZE1};
  
  /// Constructor
  HHOST_DEV ExprT(DATA* data) : m_exdata(data) {}
  
  /// Destructor
  HHOST_DEV ~ExprT(){}
  
  /// Accessor to individual entry
  HHOST_DEV typename ARG::TYPE at(size_t i) const {return m_exdata->at(i);}
  
  /// Size of the array
  HHOST_DEV size_t size() const {return m_exdata->size();}
  
  /// Accessor to the local data
  HHOST_DEV DATA* getData() const {return m_exdata;}
  
private:
  
  /// local data
  DATA* m_exdata;
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
#define EET_BINARY(__OpName__,__op__) \
template <typename V1, typename V2>			\
class __OpName__ : public ExprT< __OpName__<V1,V2>, ETPLVEC(V1,V2)> {	\
public:						\
  HHOST_DEV __OpName__ (EETYPE(V1) v1, EETYPE(V2) v2) :	    \
    ExprT<__OpName__<V1,V2>, ETPLVEC(V1,V2)>(this), e1(v1), e2(v2) {}	\
    									\
    HHOST_DEV ETPLTYPE(V1) at(size_t i) const {return __op__ ;} \
    									\
    HHOST_DEV size_t size()   const {return e1.size();}		\
 private:								\
    EETYPE(V1) e1;							\
    EETYPE(V2) e2;							\
};
  
EET_BINARY(AddT, e1.at(i)+e2.at(i))
EET_BINARY(SubT, e1.at(i)-e2.at(i))
EET_BINARY(MulT, e1.at(i)*e2.at(i))
EET_BINARY(DivT, e1.at(i)/e2.at(i))
EET_BINARY(MaxT, std::max(e1.at(i) COMMA e2.at(i)))
EET_BINARY(MinT, std::min(e1.at(i) COMMA e2.at(i)))

#undef EET_BINARY
  
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
#define EET_BINARY_C(__a1__,__a2__,__e__,__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public ExprT< __OpName__<V>, ETPL(V)> { \
public:						\
  HHOST_DEV __OpName__ (__a1__ v1, __a2__ v2) :	    \
    ExprT<__OpName__<V>, ETPL(V)>(this), e1(v1), e2(v2) {}	\
    									\
    HHOST_DEV ETPLTYPE(V) at(size_t i) const {return __op__ ;} \
    									\
    HHOST_DEV size_t size()   const {return __e__ .size();}		\
 private:								\
    __a1__ e1;								\
    __a2__ e2;								\
};
  
EET_BINARY_C(EETYPE(V), ETPLTYPE(V), e1, AddTv1, e1.at(i)+e2)
EET_BINARY_C(EETYPE(V), ETPLTYPE(V), e1, SubTv1, e1.at(i)-e2)
EET_BINARY_C(EETYPE(V), ETPLTYPE(V), e1, MulTv1, e1.at(i)*e2)
EET_BINARY_C(EETYPE(V), ETPLTYPE(V), e1, DivTv1, e1.at(i)/e2)
EET_BINARY_C(EETYPE(V), ETPLTYPE(V), e1, MaxTv1, std::max(e1.at(i) COMMA e2))
EET_BINARY_C(EETYPE(V), ETPLTYPE(V), e1, MinTv1, std::min(e1.at(i) COMMA e2))

EET_BINARY_C(ETPLTYPE(V), EETYPE(V), e2, AddTv2, e1+e2.at(i))
EET_BINARY_C(ETPLTYPE(V), EETYPE(V), e2, SubTv2, e1-e2.at(i))
EET_BINARY_C(ETPLTYPE(V), EETYPE(V), e2, MulTv2, e1*e2.at(i))
EET_BINARY_C(ETPLTYPE(V), EETYPE(V), e2, DivTv2, e1/e2.at(i))
EET_BINARY_C(ETPLTYPE(V), EETYPE(V), e2, MaxTv2, std::max(e1 COMMA e2.at(i)))
EET_BINARY_C(ETPLTYPE(V), EETYPE(V), e2, MinTv2, std::min(e1 COMMA e2.at(i)))

#undef EET_BINARY_C
  
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

#define EET_BINARY_OP(__t1__,__t2__,__a1__,__a2__,__OpName__,__op__)	\
template __t1__							\
HHOST_DEV  inline __OpName__ __t2__ __op__ (__a1__ v1, __a2__ v2) \
{									\
  return __OpName__ __t2__(v1,v2);						\
}

/// binary operator overloading taking two expressions
EET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, EETYPE(V1), EETYPE(V2), AddT, operator+)
EET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, EETYPE(V1), EETYPE(V2), SubT, operator-)
EET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, EETYPE(V1), EETYPE(V2), MulT, operator*)
EET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, EETYPE(V1), EETYPE(V2), DivT, operator/)
EET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, EETYPE(V1), EETYPE(V2), MaxT, max)
EET_BINARY_OP(<typename V1 COMMA typename V2>, <V1 COMMA V2>, EETYPE(V1), EETYPE(V2), MinT, min)

/// binary operator overloading taking one expression and one constant value
EET_BINARY_OP(<typename V>, <V>, EETYPE(V), ETPLTYPE(V), AddTv1, operator+)
EET_BINARY_OP(<typename V>, <V>, EETYPE(V), ETPLTYPE(V), SubTv1, operator-)
EET_BINARY_OP(<typename V>, <V>, EETYPE(V), ETPLTYPE(V), MulTv1, operator*)
EET_BINARY_OP(<typename V>, <V>, EETYPE(V), ETPLTYPE(V), DivTv1, operator/)
EET_BINARY_OP(<typename V>, <V>, EETYPE(V), ETPLTYPE(V), MaxTv1, max)
EET_BINARY_OP(<typename V>, <V>, EETYPE(V), ETPLTYPE(V), MinTv1, min)
  
/// binary operator overloading taking one constant value and one expression
EET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), EETYPE(V), AddTv2, operator+)
EET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), EETYPE(V), SubTv2, operator-)
EET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), EETYPE(V), MulTv2, operator*)
EET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), EETYPE(V), DivTv2, operator/)
EET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), EETYPE(V), MaxTv2, max)
EET_BINARY_OP(<typename V>, <V>, ETPLTYPE(V), EETYPE(V), MinTv2, min)
  
#undef EET_BINARY_OP 

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
#define EET_UNARY(__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public ExprT< __OpName__<V>, ETPL(V)> { \
public:						\
  HHOST_DEV __OpName__ (EETYPE(V) v1) :	    \
    ExprT<__OpName__<V>, ETPL(V)>(this), e1(v1) {}	\
    									\
    HHOST_DEV ETPLTYPE(V) at(size_t i) const {return __op__(e1.at(i));} \
    									\
    HHOST_DEV size_t size() const {return e1.size();}			\
 private:								\
    EETYPE(V) e1;							\
};

  EET_UNARY(CosT,     std::cos)
  EET_UNARY(SinT,     std::sin)
  EET_UNARY(TanT,     std::tan)
  EET_UNARY(AcosT,    std::acos)
  EET_UNARY(AsinT,    std::asin)
  EET_UNARY(AtanT,    std::atan)
  EET_UNARY(CoshT,    std::cosh)
  EET_UNARY(SinhT,    std::sinh)
  EET_UNARY(TanhT,    std::tanh)
  EET_UNARY(LogT,     std::log)
  EET_UNARY(Log10T,   std::log10)
  EET_UNARY(ExpT,     std::exp)
  EET_UNARY(SqrtT,    std::sqrt)
  EET_UNARY(AbsT,     std::abs)
  EET_UNARY(MinusOpT, -)
#undef EET_UNARY
  
//////////////////////////////////////////////////////////////////////////////
  
#define EET_UNARY_OP(__OpName__,__op__)		\
template <typename V>						\
HHOST_DEV inline __OpName__<V> __op__ (EETYPE(V) v1)		\
{									\
  return __OpName__<V>(v1);						\
}
    
  /// unary operator overloading taking one expression
  /// (handle will care: they occasionally might not work on GPU)
  EET_UNARY_OP(SinT,     sin)
  EET_UNARY_OP(CosT,     cos)
  EET_UNARY_OP(TanT,     tan) 
  EET_UNARY_OP(AsinT,    asin)
  EET_UNARY_OP(AcosT,    acos)
  EET_UNARY_OP(AtanT,    atan) 
  EET_UNARY_OP(SinhT,    sinh)
  EET_UNARY_OP(CoshT,    cosh)
  EET_UNARY_OP(TanhT,    tanh) 
  EET_UNARY_OP(LogT,     log)
  EET_UNARY_OP(Log10T,   log10) 
  EET_UNARY_OP(ExpT,     exp)
  EET_UNARY_OP(SqrtT,    sqrt) 
  EET_UNARY_OP(AbsT,     abs)
  EET_UNARY_OP(MinusOpT, operator-)
#undef EET_UNARY_OP 

//////////////////////////////////////////////////////////////////////////////
 
  } // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
