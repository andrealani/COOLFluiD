#ifndef COOLFluiD_FluctSplit_CIntegralTriagP2O3_hh
#define COOLFluiD_FluctSplit_CIntegralTriagP2O3_hh

#include "FluctSplit/MetaSchemes/ContourIntegral.hh"
#include "FluctSplit/MetaSchemes/TriagP2.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

/// Specialization for CFPolyOrder::ORDER3 contour integration of a P2 triangle
template <> struct CIntegral < TriagP2, CFPolyOrder::ORDER3 >
{
  typedef TriagP2 Elem;
  typedef CIntegral < TriagP2, CFPolyOrder::ORDER3 > Integrator;

  enum { NBSF       = Elem::NBNODES };
  enum { NBQDPT     = 27 };
  enum { NBSUBQDPT  = 9 };
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
  };

  /// weights should sum to 1.0
  static inline CFreal wl () { return 2.5 / 9.; }
  static inline CFreal wc () { return 4.0 / 9.; }
  static inline CFreal wr () { return 2.5 / 9.; }

  static inline CFreal pl () { return ( 1.0 - std::sqrt(0.6) ) * 0.250; }
  static inline CFreal pc () { return 0.250; }
  static inline CFreal pr () { return ( 1.0 + std::sqrt(0.6) ) * 0.250; }

  static inline CFreal a () { return 0.00; }
  static inline CFreal b () { return pl(); }
  static inline CFreal c () { return pc(); }
  static inline CFreal d () { return pr(); }
  static inline CFreal e () { return 0.50; }
  static inline CFreal f () { return 0.50 + pl(); }
  static inline CFreal g () { return 0.50 + pc(); }
  static inline CFreal h () { return 0.50 + pr(); }
  static inline CFreal i () { return 1.00; }

  template < unsigned int SUBELEM, unsigned int NODAL_NORMAL, typename FSDATA >
  static const RealVector& SubFaceNormal ( const FSDATA& data )
  { return data.sub_elem[SUBELEM].nodal_normals[NODAL_NORMAL]; }

  template < unsigned int SUBELEM, unsigned int NODAL_NORMAL, typename FSDATA >
  static double SubFaceLenght ( const FSDATA& data )
  { return data.sub_elem[SUBELEM].nodal_areas[NODAL_NORMAL]; }

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
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<0>::weight  () { return wl(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<0>::xi      () { return b(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<0>::eta     () { return a(); }

// definition of element quadrature point 1
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<1>::weight  () { return wc(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<1>::xi      () { return c(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<1>::eta     () { return a(); }

// definition of element quadrature point 2
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<2>::weight  () { return wr(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<2>::xi      () { return d(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<2>::eta     () { return a(); }

// definition of element quadrature point 3
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<3>::weight  () { return wl(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<3>::xi      () { return f(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<3>::eta     () { return a(); }

// definition of element quadrature point 4
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<4>::weight  () { return wc(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<4>::xi      () { return g(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<4>::eta     () { return a(); }

// definition of element quadrature point 5
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<5>::weight  () { return wr(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<5>::xi      () { return h(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<5>::eta     () { return a(); }


///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 6
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<6>::weight  () { return wl(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<6>::xi      () { return a(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<6>::eta     () { return b(); }

// definition of element quadrature point 7
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<7>::weight  () { return wl() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<7>::xi      () { return d(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<7>::eta     () { return b(); }

// definition of element quadrature point 8
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<8>::weight  () { return wl(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<8>::xi      () { return e(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<8>::eta     () { return b(); }

// definition of element quadrature point 9
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<9>::weight  () { return wl() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<9>::xi      () { return h(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<9>::eta     () { return b(); }

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 10
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<10>::weight  () { return wc(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<10>::xi      () { return a(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<10>::eta     () { return c(); }

// definition of element quadrature point 11
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<11>::weight  () { return wc() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<11>::xi      () { return c(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<11>::eta     () { return c(); }

// definition of element quadrature point 12
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<12>::weight  () { return wc(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<12>::xi      () { return e(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<12>::eta     () { return c(); }

// definition of element quadrature point 13
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<13>::weight  () { return wc() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<13>::xi      () { return g(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<13>::eta     () { return c(); }

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 14
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<14>::weight  () { return wr(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<14>::xi      () { return a(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<14>::eta     () { return d(); }

// definition of element quadrature point 15
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<15>::weight  () { return wr() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<15>::xi      () { return b(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<15>::eta     () { return d(); }

// definition of element quadrature point 16
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<16>::weight  () { return wr(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<16>::xi      () { return e(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<16>::eta     () { return d(); }

// definition of element quadrature point 17
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<17>::weight  () { return wr() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<17>::xi      () { return f(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<17>::eta     () { return d(); }

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 18
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<18>::weight  () { return wl(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<18>::xi      () { return b(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<18>::eta     () { return e(); }

// definition of element quadrature point 19
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<19>::weight  () { return wc(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<19>::xi      () { return c(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<19>::eta     () { return e(); }

// definition of element quadrature point 20
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<20>::weight  () { return wr(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<20>::xi      () { return d(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<20>::eta     () { return e(); }

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 21
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<21>::weight  () { return wl(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<21>::xi      () { return a(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<21>::eta     () { return f(); }

// definition of element quadrature point 22
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<22>::weight  () { return wl() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<22>::xi      () { return d(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<22>::eta     () { return f(); }

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 23
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<23>::weight  () { return wc(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<23>::xi      () { return a(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<23>::eta     () { return g(); }

// definition of element quadrature point 24
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<24>::weight  () { return wc() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<24>::xi      () { return c(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<24>::eta     () { return g(); }

///////////////////////////////////////////////////////////////////////////////

// definition of element quadrature point 25
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<25>::weight  () { return wr(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<25>::xi      () { return a(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<25>::eta     () { return h(); }

// definition of element quadrature point 26
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<26>::weight  () { return wr() * std::sqrt(2.0); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<26>::xi      () { return b(); }
template <> double CIntegral<TriagP2,CFPolyOrder::ORDER3>::QdPt<26>::eta     () { return h(); }

///////////////////////////////////////////////////////////////////////////////

/// Definition of subelem 0 quadrature point 0 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<0>
{
  enum { ID = 0 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 1 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<1>
{
  enum { ID = 1 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 2 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<2>
{
  enum { ID = 2 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 3 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<3>
{
  enum { ID = 7 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 4 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<4>
{
  enum { ID = 11 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 5 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<5>
{
  enum { ID = 15 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 6 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<6>
{
  enum { ID = 14 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 7 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<7>
{
  enum { ID = 10 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 0 quadrature point 8 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<0>::QdPt<8>
{
  enum { ID = 6 };
  enum { SUBELEM = 0 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

///////////////////////////////////////////////////////////////////////////////

/// Definition of subelem 1 quadrature point 0 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<0>
{
  enum { ID = 3 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 1 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<1>
{
  enum { ID = 4 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 2 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<2>
{
  enum { ID = 5 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 3 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<3>
{
  enum { ID = 9 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 4 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<4>
{
  enum { ID = 13 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 5 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<5>
{
  enum { ID = 17 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 6 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<6>
{
  enum { ID = 16 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 7 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<7>
{
  enum { ID = 12 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 1 quadrature point 8 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<1>::QdPt<8>
{
  enum { ID = 8 };
  enum { SUBELEM = 1 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

///////////////////////////////////////////////////////////////////////////////

/// Definition of subelem 2 quadrature point 0 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<0>
{
  enum { ID = 18 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 1 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<1>
{
  enum { ID = 19 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 2 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<2>
{
  enum { ID = 20 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 3 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<3>
{
  enum { ID = 22 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 4 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<4>
{
  enum { ID = 24 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 5 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<5>
{
  enum { ID = 26 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 6 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<6>
{
  enum { ID = 25 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 7 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<7>
{
  enum { ID = 23 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 2 quadrature point 8 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<2>::QdPt<8>
{
  enum { ID = 21 };
  enum { SUBELEM = 2 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

///////////////////////////////////////////////////////////////////////////////

/// Definition of subelem 3 quadrature point 0 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<0>
{
  enum { ID = 15 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) * SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 1 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<1>
{
  enum { ID = 11 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) *SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 2 (subface 0)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<2>
{
  enum { ID = 7 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 2 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return  SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return (1./std::sqrt(2.0)) *SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 3 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<3>
{
  enum { ID = 8 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 4 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<4>
{
  enum { ID = 12 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 5 (subface 1)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<5>
{
  enum { ID = 16 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 0 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 6 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<6>
{
  enum { ID = 20 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return  SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 7 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<7>
{
  enum { ID = 19 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

/// Definition of subelem 3 quadrature point 8 (subface 2)
template <> template <>
struct CIntegral<TriagP2,CFPolyOrder::ORDER3>::SubElem<3>::QdPt<8>
{
  enum { ID = 18 };
  enum { SUBELEM = 3 };
  enum { NODAL_NORMAL = 1 };
  template < typename FSDATA >
  static const RealVector& normal ( const FSDATA& data ) { return SubFaceNormal<SUBELEM,NODAL_NORMAL>(data); }
  template < typename FSDATA >
  static double face_jacob ( const FSDATA& data )  { return SubFaceLenght<SUBELEM,NODAL_NORMAL>(data); }
};

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_CIntegralTriagP2O3_hh
