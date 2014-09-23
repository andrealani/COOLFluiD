// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_ExprT_hh
#define COOLFluiD_CudaEnv_ExprT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CUDA/MacrosET.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
  
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
  HOST_DEVICE ExprT(DATA* data) : m_exdata(data) {}
  
  /// Destructor
  HOST_DEVICE ~ExprT(){}
  
  /// Accessor to individual entry
  HOST_DEVICE typename ARG::TYPE at(size_t i) const {return m_exdata->at(i);}
  
  /// Size of the array
  HOST_DEVICE size_t size() const {return m_exdata->size();}
  
  /// Accessor to the local data
  HOST_DEVICE DATA* getData() const {return m_exdata;}
  
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
#define ET_BINARY(__OpName__,__op__) \
template <typename V1, typename V2>			\
class __OpName__ : public ExprT< __OpName__<V1,V2>, TPLVEC(V1,V2)> {	\
public:						\
  HOST_DEVICE __OpName__ (ETYPE(V1) v1, ETYPE(V2) v2) :	    \
    ExprT<__OpName__<V1,V2>, TPLVEC(V1,V2)>(this), e1(v1), e2(v2) {}	\
    									\
    HOST_DEVICE TPLTYPE(V1) at(size_t i) const {return __op__ ;} \
    									\
    HOST_DEVICE size_t size()   const {return e1.size();}		\
 private:								\
    ETYPE(V1) e1;							\
    ETYPE(V2) e2;							\
};
  
ET_BINARY(AddT, e1.at(i)+e2.at(i))
ET_BINARY(SubT, e1.at(i)-e2.at(i))
ET_BINARY(MulT, e1.at(i)*e2.at(i))
ET_BINARY(DivT, e1.at(i)/e2.at(i))
ET_BINARY(MaxT, max(e1.at(i) COMA e2.at(i)))
ET_BINARY(MinT, min(e1.at(i) COMA e2.at(i)))

#undef ET_BINARY
  
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
#define ET_BINARY_C(__a1__,__a2__,__e__,__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public ExprT< __OpName__<V>, TPL(V)> { \
public:						\
  HOST_DEVICE __OpName__ (__a1__ v1, __a2__ v2) :	    \
    ExprT<__OpName__<V>, TPL(V)>(this), e1(v1), e2(v2) {}	\
    									\
    HOST_DEVICE TPLTYPE(V) at(size_t i) const {return __op__ ;} \
    									\
    HOST_DEVICE size_t size()   const {return __e__ .size();}		\
 private:								\
    __a1__ e1;								\
    __a2__ e2;								\
};
  
ET_BINARY_C(ETYPE(V), TPLTYPE(V), e1, AddTv1, e1.at(i)+e2)
ET_BINARY_C(ETYPE(V), TPLTYPE(V), e1, SubTv1, e1.at(i)-e2)
ET_BINARY_C(ETYPE(V), TPLTYPE(V), e1, MulTv1, e1.at(i)*e2)
ET_BINARY_C(ETYPE(V), TPLTYPE(V), e1, DivTv1, e1.at(i)/e2)
ET_BINARY_C(ETYPE(V), TPLTYPE(V), e1, MaxTv1, max(e1.at(i) COMA e2))
ET_BINARY_C(ETYPE(V), TPLTYPE(V), e1, MinTv1, min(e1.at(i) COMA e2))

ET_BINARY_C(TPLTYPE(V), ETYPE(V), e2, AddTv2, e1+e2.at(i))
ET_BINARY_C(TPLTYPE(V), ETYPE(V), e2, SubTv2, e1-e2.at(i))
ET_BINARY_C(TPLTYPE(V), ETYPE(V), e2, MulTv2, e1*e2.at(i))
ET_BINARY_C(TPLTYPE(V), ETYPE(V), e2, DivTv2, e1/e2.at(i))
ET_BINARY_C(TPLTYPE(V), ETYPE(V), e2, MaxTv2, max(e1 COMA e2.at(i)))
ET_BINARY_C(TPLTYPE(V), ETYPE(V), e2, MinTv2, min(e1 COMA e2.at(i)))

#undef ET_BINARY_C
  
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

#define ET_BINARY_OP(__t1__,__t2__,__a1__,__a2__,__OpName__,__op__)	\
template __t1__							\
HOST_DEVICE  inline __OpName__ __t2__ __op__ (__a1__ v1, __a2__ v2) \
{									\
  return __OpName__ __t2__(v1,v2);						\
}

/// binary operator overloading taking two expressions
ET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, ETYPE(V1), ETYPE(V2), AddT, operator+)
ET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, ETYPE(V1), ETYPE(V2), SubT, operator-)
ET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, ETYPE(V1), ETYPE(V2), MulT, operator*)
ET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, ETYPE(V1), ETYPE(V2), DivT, operator/)
ET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, ETYPE(V1), ETYPE(V2), MaxT, max)
ET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, ETYPE(V1), ETYPE(V2), MinT, min)

/// binary operator overloading taking one expression and one constant value
ET_BINARY_OP(<typename V>, <V>, ETYPE(V), TPLTYPE(V), AddTv1, operator+)
ET_BINARY_OP(<typename V>, <V>, ETYPE(V), TPLTYPE(V), SubTv1, operator-)
ET_BINARY_OP(<typename V>, <V>, ETYPE(V), TPLTYPE(V), MulTv1, operator*)
ET_BINARY_OP(<typename V>, <V>, ETYPE(V), TPLTYPE(V), DivTv1, operator/)
ET_BINARY_OP(<typename V>, <V>, ETYPE(V), TPLTYPE(V), MaxTv1, max)
ET_BINARY_OP(<typename V>, <V>, ETYPE(V), TPLTYPE(V), MinTv1, min)
  
/// binary operator overloading taking one constant value and one expression
ET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), ETYPE(V), AddTv2, operator+)
ET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), ETYPE(V), SubTv2, operator-)
ET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), ETYPE(V), MulTv2, operator*)
ET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), ETYPE(V), DivTv2, operator/)
ET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), ETYPE(V), MaxTv2, max)
ET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), ETYPE(V), MinTv2, min)
  
#undef ET_BINARY_OP 

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
#define ET_UNARY(__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public ExprT< __OpName__<V>, TPL(V)> { \
public:						\
  HOST_DEVICE __OpName__ (ETYPE(V) v1) :	    \
    ExprT<__OpName__<V>, TPL(V)>(this), e1(v1) {}	\
    									\
    HOST_DEVICE TPLTYPE(V) at(size_t i) const {return __op__(e1.at(i));} \
    									\
    HOST_DEVICE size_t size() const {return e1.size();}			\
 private:								\
    ETYPE(V) e1;							\
};

  ET_UNARY(CosT,     cos)
  ET_UNARY(SinT,     sin)
  ET_UNARY(TanT,     tan)
  ET_UNARY(AcosT,    acos)
  ET_UNARY(AsinT,    asin)
  ET_UNARY(AtanT,    atan)
  ET_UNARY(CoshT,    cosh)
  ET_UNARY(SinhT,    sinh)
  ET_UNARY(TanhT,    tanh)
  ET_UNARY(LogT,     log)
  ET_UNARY(Log10T,   log10)
  ET_UNARY(ExpT,     exp)
  ET_UNARY(SqrtT,    sqrt)
  ET_UNARY(AbsT,     abs)
  ET_UNARY(MinusOpT, -)
#undef ET_UNARY
  
//////////////////////////////////////////////////////////////////////////////
  
#define ET_UNARY_OP(__OpName__,__op__)		\
template <typename V>						\
HOST_DEVICE inline __OpName__<V> __op__ (ETYPE(V) v1)		\
{									\
  return __OpName__<V>(v1);						\
}
    
  /// unary operator overloading taking one expression
  /// (handle will care: they occasionally might not work on GPU)
  ET_UNARY_OP(SinT,     Sin)
  ET_UNARY_OP(CosT,     Cos)
  ET_UNARY_OP(TanT,     Tan) 
  ET_UNARY_OP(AsinT,    Asin)
  ET_UNARY_OP(AcosT,    Acos)
  ET_UNARY_OP(AtanT,    Atan) 
  ET_UNARY_OP(SinhT,    Sinh)
  ET_UNARY_OP(CoshT,    Cosh)
  ET_UNARY_OP(TanhT,    Tanh) 
  ET_UNARY_OP(LogT,     Log)
  ET_UNARY_OP(Log10T,   Log10) 
  ET_UNARY_OP(ExpT,     Exp)
  ET_UNARY_OP(SqrtT,    Sqrt) 
  ET_UNARY_OP(AbsT,     Abs)
  ET_UNARY_OP(MinusOpT, operator-)
#undef ET_UNARY_OP 

//////////////////////////////////////////////////////////////////////////////
 
  } // namespace CudaEnv

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
