// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MacrosET_hh
#define COOLFluiD_MathTools_MacrosET_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

#define COMMA ,

//////////////////////////////////////////////////////////////////////////////

#define EEREF(a) const a&
#define HHOST_DEV

//////////////////////////////////////////////////////////////////////////////

#define EVECLOOP(idx,imin,imax,arg) for (size_t idx = imin; idx < imax; ++idx) arg
#define EMATLOOP(i,imin,imax,j, jmin,jmax,arg)				\
  for (size_t i = imin; i < imax; ++i) for (size_t j = jmin; j < jmax; ++j) arg

//////////////////////////////////////////////////////////////////////////////

#define ENMAX(a,b) ECMP<a,b>::MAX
#define EGETSIZE1(a) ((a > 0) ? a : size())
#define EGETSIZER(a) ((a > 0) ? a : nbRows())
#define EGETSIZEC(a) ((a > 0) ? a : nbCols())
#define EGETSIZE12(a,b) ((a > 0 && b > 0) ? a*b : nbRows()*nbCols())

//////////////////////////////////////////////////////////////////////////////

// Macros for tuple manipulation 
#define ETPT(a) a::TUPLE
#define ETPL(a) typename a::TUPLE
#define ETPLTYPE(a) typename a::TUPLE::TYPE
#define ETPLVEC(a,b)     ETuple<ETPLTYPE(a),ENMAX(a::SIZE1, b::SIZE1),0>
#define ETPLMAT(a,b)     ETuple<ETPLTYPE(a),ENMAX(a::SIZE1, b::SIZE1),ENMAX(a::SIZE2, b::SIZE2)>
#define ETPLMATMULT(a,b) ETuple<ETPLTYPE(a),ENMAX(a::SIZE1, 0),ENMAX(b::SIZE2, 0)>

//////////////////////////////////////////////////////////////////////////////

#define EETYPE(a)   EEREF(ExprT<a COMMA ETPL(a)>)
#define EETYPEV(a)  EEREF(ExprT<a COMMA typename a::TUPLE>)
#define MEETYPE(a)  EEREF(MatExprT<a COMMA ETPL(a)>)
#define MEETYPEV(a) EEREF(MatExprT<a COMMA typename a::TUPLE>)
  
//////////////////////////////////////////////////////////////////////////////

/// This struct determines at compile time the bigger integer between N and M
/// @author Andrea Lani
template <int N, int M> struct ECMP {enum {MAX=(N>M) ? N : M};};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

}   // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
