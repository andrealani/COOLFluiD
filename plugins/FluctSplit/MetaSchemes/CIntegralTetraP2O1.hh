#ifndef COOLFluiD_FluctSplit_CIntegralTetraP2O1_hh
#define COOLFluiD_FluctSplit_CIntegralTetraP2O1_hh

#include "FluctSplit/MetaSchemes/ContourIntegral.hh"
#include "FluctSplit/MetaSchemes/TetraP2.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// Specialization for first order contour integration of a triangle
template <> struct CIntegral < TetraP2, CFPolyOrder::ORDER1 >
{
  typedef TetraP2 Elem;
  typedef CIntegral < TetraP2, CFPolyOrder::ORDER1 > Integrator;

  enum { NBSF       = Elem::NBNODES };
  enum { NBQDPT     = 4 };
  enum { NBSUBQDPT  = 4 };
  enum { NBSUBELEM  = Elem::NBSUBELEM  };
  enum { NBSUBNODES = Elem::NBSUBNODES };

  /// Class for each quadrature point in the element
  /// This template must be specialized.
  template < unsigned int QDPT >
  struct QdPt
  {
    /// @post Sum of the weights of all quad points in a subelement face is 1.
    static double weight ();
    static double xi     ();
    static double eta    ();
    static double zeta    ();
  };

  /// Class for each quadrature point
  /// This template must be specialized.
  template < unsigned int NELEM >
  struct SubElem
  {
    /// Class for returning the normals as seen form the subelement quadrature point
    template < unsigned int QDPT >
    struct QdPt
    {
      /// Returns the normal on the subelem quadrature point,
      /// as seen from the subelement.
      /// @post normals points to the inside of the subelement and is a unit normal with size 1
      /// @return RealVector of size DIM with the normal components
      template < typename FSDATA >
      static const RealVector& normal ( const FSDATA& data );
      /// Returns the face jacobian at the quadrature point
      /// @post On a straight face it is the ratio between the face area
      ///       and the reference face area
      /// @return RealVector of size DIM with the normal components
      template < typename FSDATA >
      static double face_jacob ( const FSDATA& data );
    };
  };
};

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 0
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<0>::weight  () { return 0.500; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<0>::xi      () { return 1./3.; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<0>::eta     () { return 1./3.; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<0>::zeta    () { return 0.000; }

// definition of element quadrature point 1
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<1>::weight  () { return 0.500; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<1>::xi      () { return 1./3.; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<1>::eta     () { return 0.000; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<1>::zeta    () { return 1./3.; }

// definition of element quadrature point 2
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<2>::weight  () { return std::sqrt(3.0)/2.0; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<2>::xi      () { return 1./3.; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<2>::eta     () { return 1./3.; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<2>::zeta    () { return 1./3.; }

// definition of element quadrature point 3
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<3>::weight  () { return 0.500; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<3>::xi      () { return 0.000; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<3>::eta     () { return 1./3.; }
template <> double CIntegral<TetraP2,CFPolyOrder::ORDER1>::QdPt<3>::zeta    () { return 1./3.; }

///////////////////////////////////////////////////////////////////////////////

/// Definition of subelem 0 quadrature point 0
template <> template <>
struct CIntegral<TetraP2,CFPolyOrder::ORDER1>::SubElem<0>::QdPt<0>
{
  enum { ID = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return data.sub_elem[0].nodal_normals[3]; }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return data.sub_elem[0].nodal_areas[3]; }
};

/// Definition of subelem 0 quadrature point 1
template <> template <>
struct CIntegral<TetraP2,CFPolyOrder::ORDER1>::SubElem<0>::QdPt<1>
{
  enum { ID = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return data.sub_elem[0].nodal_normals[2]; }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return data.sub_elem[0].nodal_areas[2]; }
};

/// Definition of subelem 0 quadrature point 2
template <> template <>
struct CIntegral<TetraP2,CFPolyOrder::ORDER1>::SubElem<0>::QdPt<2>
{
  enum { ID = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return data.sub_elem[0].nodal_normals[0]; }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return ( std::sqrt(3.0) / 2.0 ) * data.sub_elem[0].nodal_areas[0]; }
};

/// Definition of subelem 0 quadrature point 3
template <> template <>
struct CIntegral<TetraP2,CFPolyOrder::ORDER1>::SubElem<0>::QdPt<3>
{
  enum { ID = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return data.sub_elem[0].nodal_normals[1]; }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return data.sub_elem[0].nodal_areas[1]; }
};

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_CIntegralTetraP2O1_hh
