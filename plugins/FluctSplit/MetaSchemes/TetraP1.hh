#ifndef COOLFluiD_FluctSplit_TetraP1_hh
#define COOLFluiD_FluctSplit_TetraP1_hh

#include "FluctSplit/MetaSchemes/BaseSubElem.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

struct TetraP1
{
  enum { SHAPE      = CFGeoShape::TETRA };
  enum { DIM        = DIM_3D };
  enum { ORDER      = CFPolyOrder::ORDER1 };
  enum { NBNODES    = 4 };
  enum { NBSUBELEM  = 1 };
  enum { NBSUBNODES = 4 };

  static std::string getClassName() { return "TetraP1"; };
  static std::string getShapeName() { return "Tetra"; };

  /// Class for each shape function, template must be specialized.
  template < unsigned int NODE >
  struct ShapeF
  {
    /// Gives value of each shape funtion, the coordinates are a template parameter,
    /// that should be a class providing the functions xi(), eta() and zeta if in 3D.
    /// @return value of the function at the coordinates
    template < typename COORD > static double value ();
  };

  template < unsigned int NELEM >
  struct SubElem : public BaseSubElem<DIM>
  {
    /// Template class that maps id of subelem node to node in master element
    template < unsigned int NODE > struct SNode {};
    /// Computes the normals in the subelement and places then in the FSDATA
    template < typename FSDATA >
    static void compute_nodal_normals ( FSDATA& data);
  };

};

///////////////////////////////////////////////////////////////////////////////

// Nodes of subelement 0
template <> template <> struct TetraP1::SubElem<0>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TetraP1::SubElem<0>::SNode<1> { enum { ID = 1 }; };
template <> template <> struct TetraP1::SubElem<0>::SNode<2> { enum { ID = 2 }; };
template <> template <> struct TetraP1::SubElem<0>::SNode<3> { enum { ID = 3 }; };

///////////////////////////////////////////////////////////////////////////////

/// shape function node 0 \f[ N = 1 - \ksi - \eta ) \f]
template <> template < typename COORD >
double TetraP1::ShapeF<0>::value () { return 1.0 - ( COORD::xi() + COORD::eta() + COORD::zeta()); }

/// shape function node 1 \f[ N = \xi \f]
template <> template < typename COORD >
double TetraP1::ShapeF<1>::value () { return COORD::xi(); }

/// shape function node 2 \f[ N = \eta \f]
template <> template < typename COORD >
double TetraP1::ShapeF<2>::value () { return COORD::eta(); }

/// shape function node 3 \f[ N = \zeta \f]
template <> template < typename COORD >
double TetraP1::ShapeF<3>::value () { return COORD::zeta(); }

///////////////////////////////////////////////////////////////////////////////

template < unsigned int NELEM >
template < typename FSDATA >
void TetraP1::SubElem<NELEM>::compute_nodal_normals ( FSDATA& data )
{
  std::vector<RealVector>& nodal_normals = data.sub_elem[NELEM].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[NELEM].nodal_areas;

  typedef typename TetraP1::SubElem<NELEM>::template SNode<0> SNode0;
  typedef typename TetraP1::SubElem<NELEM>::template SNode<1> SNode1;
  typedef typename TetraP1::SubElem<NELEM>::template SNode<2> SNode2;
  typedef typename TetraP1::SubElem<NELEM>::template SNode<3> SNode3;

  compute_normal_vector( *data.nodes[SNode0::ID], *data.nodes[SNode1::ID], *data.nodes[SNode2::ID], nodal_normals[0]); // Face 021 inverted
  compute_normal_vector( *data.nodes[SNode0::ID], *data.nodes[SNode3::ID], *data.nodes[SNode1::ID], nodal_normals[1]); // Face 013 inverted
  compute_normal_vector( *data.nodes[SNode1::ID], *data.nodes[SNode3::ID], *data.nodes[SNode2::ID], nodal_normals[2]); // Face 123 inverted
  compute_normal_vector( *data.nodes[SNode0::ID], *data.nodes[SNode2::ID], *data.nodes[SNode3::ID], nodal_normals[3]); // Face 032 inverted

  // compute nodal areas
  nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
  nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
  nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
  nodal_areas[3] = nodal_normals[3].norm2(); nodal_normals[3] /= nodal_areas[3];
}

///////////////////////////////////////////////////////////////////////////////

// template <>
// template < typename FSDATA >
// void TetraP1::SubElem<0>::compute_nodal_normals ( FSDATA& data )
// {
//   std::vector<RealVector>& nodal_normals = data.sub_elem[0].nodal_normals;
//   std::vector<CFreal>&     nodal_areas   = data.sub_elem[0].nodal_areas;
//
//   // Face 021 inverted
//   compute_normal_vector( *data.nodes[0], *data.nodes[1], *data.nodes[2], nodal_normals[0]);
//
//   // Face 013 inverted
//   compute_normal_vector( *data.nodes[0], *data.nodes[3], *data.nodes[1], nodal_normals[1]);
//
//   // Face 123 inverted
//   compute_normal_vector( *data.nodes[1], *data.nodes[3], *data.nodes[2], nodal_normals[2]);
//
//   // Face 032 inverted
//   compute_normal_vector( *data.nodes[0], *data.nodes[2], *data.nodes[3], nodal_normals[3]);
//
//   // compute nodal areas
//   nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
//   nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
//   nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
//   nodal_areas[3] = nodal_normals[3].norm2(); nodal_normals[3] /= nodal_areas[3];
// }

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_TetraP1_hh
