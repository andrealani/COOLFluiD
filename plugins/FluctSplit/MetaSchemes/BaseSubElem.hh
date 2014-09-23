#ifndef COOLFluiD_FluctSplit_BaseSubElem_hh
#define COOLFluiD_FluctSplit_BaseSubElem_hh

// #include "Common/CFLog.hh"

#include "Common/StringOps.hh"
#include "Common/Meta/Power.hh"
#include "MathTools/RealVector.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/CFPolyOrder.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

template < unsigned int DIM >
struct BaseSubElem
{
  inline static void compute_normal_vector (const RealVector& a, const RealVector& b, RealVector& n);
  inline static void compute_normal_vector (const RealVector& a, const RealVector& b, const RealVector& c, RealVector& n);
};

template <>
inline void BaseSubElem<DIM_2D>::compute_normal_vector (const RealVector& a, const RealVector& b, RealVector& n)
{
  n[XX] = b[YY] - a[YY];
  n[YY] = a[XX] - b[XX];
}

template <>
inline void BaseSubElem<DIM_3D>::compute_normal_vector (const RealVector& a, const RealVector& b, const RealVector& c, RealVector& n)
{
    const CFreal v1x = b[XX] - a[XX];
    const CFreal v1y = b[YY] - a[YY];
    const CFreal v1z = b[ZZ] - a[ZZ];

    const CFreal v2x = c[XX] - b[XX];
    const CFreal v2y = c[YY] - b[YY];
    const CFreal v2z = c[ZZ] - b[ZZ];

    n[XX] = 0.5 * (  v1y*v2z - v1z*v2y );
    n[YY] = 0.5 * ( -v1x*v2z + v1z*v2x );
    n[ZZ] = 0.5 * (  v1x*v2y - v1y*v2x );
}

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_BaseSubElem_hh
