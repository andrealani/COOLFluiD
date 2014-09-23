#ifndef COOLFluiD_FluctSplit_ShapeFunctions_hh
#define COOLFluiD_FluctSplit_ShapeFunctions_hh

#include "FluctSplit/MetaSchemes/BaseSubElem.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

#if 0

struct TriagP2
{
  enum { SHAPE    = TRIAG };
  enum { DIM      = DIM_2D };
  enum { ORDER    = CFPolyOrder::ORDER2 };
  enum { NBNODES  = 6 }; // (n^2+n)/2
  enum { NBSUBELEM  = Common::Power<ORDER,DIM>::value };
  enum { NBSUBNODES = 3 };

  static std::string getClassName() { return "TriagP2"; };
  static std::string getShapeName() { return CFGeoShape::Convert::to_str(CFGeoShape::TRIAG); };

  template < unsigned int NELEM >
  struct SubElem
  {
    static const int table [];
  };
};

template <> const int TriagP2::SubElem<0>::table [TriagP2::NBSUBNODES] = { 0, 3, 5 };
template <> const int TriagP2::SubElem<1>::table [TriagP2::NBSUBNODES] = { 3, 1, 4 };
template <> const int TriagP2::SubElem<2>::table [TriagP2::NBSUBNODES] = { 5, 4, 2 };
template <> const int TriagP2::SubElem<3>::table [TriagP2::NBSUBNODES] = { 3, 4, 5 };

struct QuadP1
{
  enum { SHAPE    = QUAD };
  enum { DIM      = DIM_2D };
  enum { ORDER    = CFPolyOrder::ORDER1 };
  enum { NBNODES  = 4 };
  enum { NBSUBELEM  = 1 };
  enum { NBSUBNODES = 4 };

  static std::string getClassName() { return "QuadP1"; };
  static std::string getShapeName() { return "Quad"; };
};

struct QuadP2
{
  enum { SHAPE    = QUAD };
  enum { DIM      = DIM_2D };
  enum { ORDER    = CFPolyOrder::ORDER2 };
  enum { NBNODES  = 9 };
  enum { NBSUBELEM  = 4 };
  enum { NBSUBNODES = 4 };

  static std::string getClassName() { return "QuadP2"; };
  static std::string getShapeName() { return "Quad"; };
};

struct TetraP1
{
  enum { SHAPE    = TETRA };
  enum { DIM      = DIM_3D };
  enum { ORDER    = CFPolyOrder::ORDER1 };
  enum { NBNODES  = 4 };
  enum { NBSUBELEM  = Common::Power<ORDER,DIM>::value };
  enum { NBSUBNODES = 4 };

  static std::string getClassName() { return "TetraP1"; };
  static std::string getShapeName() { return "Tetra"; };
};

struct TetraP2
{
  enum { SHAPE    = TETRA };
  enum { DIM      = DIM_3D };
  enum { ORDER    = CFPolyOrder::ORDER2 };
  enum { NBNODES  = 10 };
  enum { NBSUBELEM  = Common::Power<ORDER,DIM>::value };
  enum { NBSUBNODES = 4 };

  static std::string getClassName() { return "TetraP2"; };
  static std::string getShapeName() { return "Tetra"; };
};
#endif

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_ShapeFunctions_hh
