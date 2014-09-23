// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_MacrosET_hh
#define COOLFluiD_CudaEnv_MacrosET_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Compatibility.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {

//////////////////////////////////////////////////////////////////////////////

#define COMA ,

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_HAVE_CUDA
#define EREF(a) const a&
#else
#define EREF(a) a
#endif 

//////////////////////////////////////////////////////////////////////////////

#define VECLOOP(idx,imin,imax,arg) for (size_t idx = imin; idx < imax; ++idx) arg
#define MATLOOP(i,imin,imax,j, jmin,jmax,arg)				\
  for (size_t i = imin; i < imax; ++i) for (size_t j = jmin; j < jmax; ++j) arg

//////////////////////////////////////////////////////////////////////////////

#define NMAX(a,b) CMP<a,b>::MAX
#define GETSIZE1(a) ((a > 0) ? a : size())
#define GETSIZER(a) ((a > 0) ? a : nbRows())
#define GETSIZEC(a) ((a > 0) ? a : nbCols())
#define GETSIZE12(a,b) ((a > 0 && b > 0) ? a*b : nbRows()*nbCols())

//////////////////////////////////////////////////////////////////////////////

// Macros for tuple manipulation 
#define TPT(a) a::TUPLE
#define TPL(a) typename a::TUPLE
#define TPLTYPE(a) typename a::TUPLE::TYPE
#define TPLVEC(a,b)     ETuple<TPLTYPE(a),NMAX(a::SIZE1, b::SIZE1),0>
#define TPLMAT(a,b)     ETuple<TPLTYPE(a),NMAX(a::SIZE1, b::SIZE1),NMAX(a::SIZE2, b::SIZE2)>
#define TPLMATMULT(a,b) ETuple<TPLTYPE(a),NMAX(a::SIZE1, 0),NMAX(b::SIZE2, 0)>

//////////////////////////////////////////////////////////////////////////////

#define ETYPE(a)   EREF(ExprT<a COMA TPL(a)>)
#define ETYPEV(a)  EREF(ExprT<a COMA typename a::TUPLE>)
#define METYPE(a)  EREF(MatExprT<a COMA TPL(a)>)
#define METYPEV(a) EREF(MatExprT<a COMA typename a::TUPLE>)
  
//////////////////////////////////////////////////////////////////////////////

/// This struct determines at compile time the bigger integer between N and M
/// @author Andrea Lani
template <int N, int M> struct CMP {enum {MAX=(N>M) ? N : M};};

//////////////////////////////////////////////////////////////////////////////

  } // namespace CudaEnv

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
