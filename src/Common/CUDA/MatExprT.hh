// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_MatExprT_hh
#define COOLFluiD_CudaEnv_MatExprT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CUDA/ExprT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
  
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
  HOST_DEVICE MatExprT(DATA* data) : ExprT<DATA, ARG>(data) {}
  
  /// Destructor
  HOST_DEVICE ~MatExprT(){} 
  
  /// Number of columns
  HOST_DEVICE size_t nbCols() const {return this->getData()->nbCols();}
  
  /// Number of rows
  HOST_DEVICE size_t nbRows() const {return this->getData()->nbRows();}

  /// Accessor to individual entry
  HOST_DEVICE typename ARG::TYPE at(size_t i, size_t j) const {return this->getData()->at(i,j);} 
  
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
#define MET_BINARY(__OpName__,__op__) \
template <typename V1, typename V2>			\
class __OpName__ : public MatExprT< __OpName__<V1,V2>, TPLMAT(V1,V2)> { \
public:						\
  HOST_DEVICE __OpName__ (METYPE(V1) v1, METYPE(V2) v2) :	    \
    MatExprT<__OpName__<V1,V2>, TPLMAT(V1,V2)>(this), e1(v1), e2(v2) {}	\
    									\
    HOST_DEVICE TPLTYPE(V1) at(size_t i, size_t j) const {return e1.at(i,j) __op__ e2.at(i,j);} \
    									\
    HOST_DEVICE size_t size()   const {return e1.size();}		\
    HOST_DEVICE size_t nbCols() const {return e1.nbCols();}        \
    HOST_DEVICE size_t nbRows() const {return e1.nbRows();}        \
									\
    METYPE(V1) e1;							\
    METYPE(V2) e2;							\
};
  
MET_BINARY(MAddT,+)
MET_BINARY(MSubT,-)
MET_BINARY(MDivT,/)

#undef MET_BINARY
  
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
#define MET_BINARY_C(__a1__,__a2__,__e1__,__e2__,__e__,__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public MatExprT< __OpName__<V>, TPL(V)> { \
public:						\
  HOST_DEVICE __OpName__ (__a1__ v1, __a2__ v2) :	    \
    MatExprT<__OpName__<V>, TPL(V)>(this), e1(v1), e2(v2) {}	\
    									\
    HOST_DEVICE TPLTYPE(V) at(size_t i, size_t j) const {return __e1__ __op__ __e2__;} \
    									\
    HOST_DEVICE size_t size()   const {return __e__ .size();}		\
    HOST_DEVICE size_t nbCols() const {return __e__ .nbCols();}		\
    HOST_DEVICE size_t nbRows() const {return __e__ .nbRows();}		\
									\
    __a1__ e1;								\
    __a2__ e2;								\
};
  
MET_BINARY_C(METYPE(V), TPLTYPE(V), e1.at(i,j), e2, e1, MAddTv1,+)
MET_BINARY_C(METYPE(V), TPLTYPE(V), e1.at(i,j), e2, e1, MSubTv1,-)
MET_BINARY_C(METYPE(V), TPLTYPE(V), e1.at(i,j), e2, e1, MMulTv1,*)
MET_BINARY_C(METYPE(V), TPLTYPE(V), e1.at(i,j), e2, e1, MDivTv1,/)

MET_BINARY_C(TPLTYPE(V), METYPE(V), e1, e2.at(i,j), e2, MAddTv2,+)
MET_BINARY_C(TPLTYPE(V), METYPE(V), e1, e2.at(i,j), e2, MSubTv2,-)
MET_BINARY_C(TPLTYPE(V), METYPE(V), e1, e2.at(i,j), e2, MMulTv2,*)
MET_BINARY_C(TPLTYPE(V), METYPE(V), e1, e2.at(i,j), e2, MDivTv2,/)

#undef MET_BINARY_C
  
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

#define MET_BINARY_OP(__t1__,__t2__,__a1__,__a2__,__OpName__,__op__)	\
template __t1__							\
HOST_DEVICE  inline __OpName__ __t2__ operator __op__ (__a1__ v1, __a2__ v2) \
{									\
  return __OpName__ __t2__(v1,v2);						\
}

/// binary operator overloading taking two expressions
MET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, METYPE(V1), METYPE(V2), MAddT, +)
MET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, METYPE(V1), METYPE(V2), MSubT, -)
MET_BINARY_OP(<typename V1 COMA typename V2>, <V1 COMA V2>, METYPE(V1), METYPE(V2), MDivT, /)

/// binary operator overloading taking one expression and one constant value
MET_BINARY_OP(<typename V>, <V>, METYPE(V), TPLTYPE(V), MAddTv1, +)
MET_BINARY_OP(<typename V>, <V>, METYPE(V), TPLTYPE(V), MSubTv1, -)
MET_BINARY_OP(<typename V>, <V>, METYPE(V), TPLTYPE(V), MMulTv1, *)
MET_BINARY_OP(<typename V>, <V>, METYPE(V), TPLTYPE(V), MDivTv1, /)
  
/// binary operator overloading taking one constant value and one expression
MET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), METYPE(V), MAddTv2, +)
MET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), METYPE(V), MSubTv2, -)
MET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), METYPE(V), MMulTv2, *)
MET_BINARY_OP(<typename V>, <V>, TPLTYPE(V), METYPE(V), MDivTv2, /)
  
#undef MET_BINARY_OP 

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
#define MET_UNARY(__OpName__,__op__) \
template <typename V>			\
class __OpName__ : public MatExprT< __OpName__<V>, TPL(V)> { \
public:						\
  HOST_DEVICE __OpName__ (METYPE(V) v1) :	    \
    MatExprT<__OpName__<V>, TPL(V)>(this), e1(v1) {}	\
    									\
    HOST_DEVICE TPLTYPE(V) at(size_t i, size_t j) const {return __op__(e1.at(i,j));} \
    									\
    HOST_DEVICE size_t size() const {return e1.size();}			\
    HOST_DEVICE size_t nbCols() const {return e1.nbCols();}        \
    HOST_DEVICE size_t nbRows() const {return e1.nbRows();}		\
    									\
    METYPE(V) e1;							\
};

  MET_UNARY(MCosT,     cos)
  MET_UNARY(MSinT,     sin)
  MET_UNARY(MTanT,     tan)
  MET_UNARY(MAcosT,    acos)
  MET_UNARY(MAsinT,    asin)
  MET_UNARY(MAtanT,    atan)
  MET_UNARY(MCoshT,    cosh)
  MET_UNARY(MSinhT,    sinh)
  MET_UNARY(MTanhT,    tanh)
  MET_UNARY(MLogT,     log)
  MET_UNARY(MLog10T,   log10)
  MET_UNARY(MExpT,     exp)
  MET_UNARY(MSqrtT,    sqrt)
  MET_UNARY(MAbsT,     abs)
  MET_UNARY(MMinusOpT, -)

#undef MET_UNARY
  
//////////////////////////////////////////////////////////////////////////////
  
#define MET_UNARY_OP(__OpName__,__op__)		\
template <typename V>						\
HOST_DEVICE inline __OpName__<V> __op__ (METYPE(V) v1)		\
{									\
  return __OpName__<V>(v1);						\
}
    
  /// unary operator overloading taking one expression
  /// (handle will care: they occasionally might not work on GPU)
  MET_UNARY_OP(MSinT,     Sin)
  MET_UNARY_OP(MCosT,     Cos)
  MET_UNARY_OP(MTanT,     Tan) 
  MET_UNARY_OP(MAsinT,    Asin)
  MET_UNARY_OP(MAcosT,    Acos)
  MET_UNARY_OP(MAtanT,    Atan) 
  MET_UNARY_OP(MSinhT,    Sinh)
  MET_UNARY_OP(MCoshT,    Cosh)
  MET_UNARY_OP(MTanhT,    Tanh) 
  MET_UNARY_OP(MLogT,     Log)
  MET_UNARY_OP(MLog10T,   Log10) 
  MET_UNARY_OP(MExpT,     Exp)
  MET_UNARY_OP(MSqrtT,    Sqrt) 
  MET_UNARY_OP(MAbsT,     Abs)
  MET_UNARY_OP(MMinusOpT, operator-)

#undef MET_UNARY_OP 

//////////////////////////////////////////////////////////////////////////////
 
/**
 * Expression class for multiplication between two CFMat-based expressions 
 * (V1 nd V2). The constructor accepts expressions V1 and V2 as arguments.
 *
 * @author Andrea Lani
 */  
template <typename V1, typename V2>
class MMulT :public MatExprT< MMulT<V1,V2>, TPLMATMULT(V1,V2)> {
public:	
  
  /// Constructor
  HOST_DEVICE MMulT (METYPE(V1) v1, METYPE(V2) v2) :
    MatExprT< MMulT<V1,V2>, TPLMATMULT(V1,V2)>(this), ex1(v1), ex2(v2) {}
    
  /// Accessor to individual entry
  HOST_DEVICE TPLTYPE(V2) at(size_t i, size_t j) const  
  { 
    TPLTYPE(V2) res = TPLTYPE(V2)();
    const size_t nbCols = ex1.nbCols();
    for (size_t k = 0; k < nbCols; ++k) {
      res += ex1.at(i,k)*ex2.at(k,j);
    }
    return res;
  }
  
  /// Size of the array
  HOST_DEVICE size_t size() const {return ex1.nbRows()*ex2.nbCols();}
  
  /// Number of columns
  HOST_DEVICE size_t nbCols() const {return ex2.nbCols();}
  
  /// Number of rows
  HOST_DEVICE size_t nbRows() const {return ex1.nbRows();}
  
private:
  
  METYPE(V1) ex1;			
  METYPE(V2) ex2;			
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename V1, typename V2>
HOST_DEVICE inline MMulT<V1,V2> operator* (METYPE(V1) v1, METYPE(V2) v2)
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
class VMMulT :public MatExprT< VMMulT<V1,V2>, TPLMATMULT(V1,V2)> {
public:	
  
  /// Constructor
  HOST_DEVICE VMMulT (ETYPE(V1) v1, METYPE(V2) v2) :
    MatExprT< VMMulT<V1,V2>, TPLMATMULT(V1,V2)>(this), ex1(v1), ex2(v2) {}
  
  /// Accessor to individual entry
  HOST_DEVICE TPLTYPE(V2) at(size_t i, size_t j) const {return ex1.at(i)*ex2.at(i,j);}
  
  /// Size of the array
  HOST_DEVICE size_t size() const {return ex2.size();}
  
  /// Number of columns
  HOST_DEVICE size_t nbCols() const {return ex2.nbCols();}
  
  /// Number of rows
  HOST_DEVICE size_t nbRows() const {return ex2.nbRows();}
 
private:
  
  ETYPE(V1)   ex1;			
  METYPE(V2)  ex2;			
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename V1, typename V2>
HOST_DEVICE inline VMMulT<V1,V2> operator* (ETYPE(V1) v1, METYPE(V2) v2)
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
class MVMulVT :public ExprT< MVMulVT<V1,V2>, TPLVEC(V1,V2)> {
public:	
  
  /// Constructor
  HOST_DEVICE MVMulVT (METYPE(V1) v1, ETYPE(V2) v2) :
    ExprT< MVMulVT<V1,V2>, TPLVEC(V1,V2)>(this), ex1(v1), ex2(v2) {}
    
  /// Accessor to individual entry
  HOST_DEVICE TPLTYPE(V2) at(size_t i) const  
  { 
    TPLTYPE(V2) res = TPLTYPE(V2)();
    for (size_t k = 0; k < ex2.size(); ++k) {
      res += ex1.at(i,k) * ex2.at(k);
    }
    return res;
  }
  
  /// Size of the array
  HOST_DEVICE size_t size() const {return ex2.size();}

private:
  
  METYPE(V1) ex1;			
  ETYPE(V2)  ex2;			
};
  
//////////////////////////////////////////////////////////////////////////////

template <typename V1, typename V2>
HOST_DEVICE inline MVMulVT<V1,V2> operator* (METYPE(V1) v1, ETYPE(V2) v2)
{
  return MVMulVT<V1,V2>(v1,v2);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace CudaEnv

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
