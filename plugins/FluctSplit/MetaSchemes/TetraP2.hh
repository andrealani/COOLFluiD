#ifndef COOLFluiD_FluctSplit_TetraP2_hh
#define COOLFluiD_FluctSplit_TetraP2_hh

#include "FluctSplit/MetaSchemes/BaseSubElem.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

struct TetraP2
{
  enum { SHAPE      = TETRA };
  enum { DIM        = DIM_3D };
  enum { ORDER      = CFPolyOrder::ORDER2 };
  enum { NBNODES    = 10 ;
  enum { NBSUBELEM  = 7 };
  enum { NBSUBNODES = 4 };

  static std::string getClassName() { return "TetraP2"; };
  static std::string getShapeName() { return "Tetra"; };

  /// Class for each shape function, template must be specialized.
  template < unsigned int NODE >
  struct ShapeF
  {
    /// computes  \f[ 1 - \ksi - \eta  - \zeta  \f]
    template < typename COORD > inline static double lambda () { return  1. -  ( COORD::xi() + COORD::eta() + COORD::zeta() ) };

    /// Gives value of each shape funtion, the coordinates are a template parameter,
    /// that should be a class providing the functions xi(), eta() and zeta if in 3D.
    /// @return value of the function at the coordinates
    template < typename COORD > inline static double value ();
  };

  template < unsigned int NELEM >
  struct SubElem : public BaseSubElem<DIM>
  {
    /// Template class that maps id of subelem node to node in master element
    template < unsigned int NODE > struct SNode {};
    /// Computes the normals in the subelement and places then in the FSDATA
    template < typename FSDATA >
    inline static void compute_nodal_normals ( FSDATA& data);
  };

};

///////////////////////////////////////////////////////////////////////////////

// Nodes of subelement 0
template <> template <> struct TetraP2::SubElem<0>::SNode<0> { enum { ID = 4 }; };
template <> template <> struct TetraP2::SubElem<0>::SNode<1> { enum { ID = 1 }; };
template <> template <> struct TetraP2::SubElem<0>::SNode<2> { enum { ID = 5 }; };
template <> template <> struct TetraP2::SubElem<0>::SNode<3> { enum { ID = 7 }; };

// Nodes of subelement 1
template <> template <> struct TetraP2::SubElem<1>::SNode<0> { enum { ID = 6 }; };
template <> template <> struct TetraP2::SubElem<1>::SNode<1> { enum { ID = 5 }; };
template <> template <> struct TetraP2::SubElem<1>::SNode<2> { enum { ID = 2 }; };
template <> template <> struct TetraP2::SubElem<1>::SNode<3> { enum { ID = 8 }; };

// Nodes of subelement 2
template <> template <> struct TetraP2::SubElem<2>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TetraP2::SubElem<2>::SNode<1> { enum { ID = 4 }; };
template <> template <> struct TetraP2::SubElem<2>::SNode<2> { enum { ID = 5 }; };
template <> template <> struct TetraP2::SubElem<2>::SNode<3> { enum { ID = 7 }; };

// Nodes of subelement 3
template <> template <> struct TetraP2::SubElem<3>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TetraP2::SubElem<3>::SNode<1> { enum { ID = 5 }; };
template <> template <> struct TetraP2::SubElem<3>::SNode<2> { enum { ID = 6 }; };
template <> template <> struct TetraP2::SubElem<3>::SNode<3> { enum { ID = 8 }; };

// Nodes of subelement 4
template <> template <> struct TetraP2::SubElem<4>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TetraP2::SubElem<4>::SNode<1> { enum { ID = 7 }; };
template <> template <> struct TetraP2::SubElem<4>::SNode<2> { enum { ID = 5 }; };
template <> template <> struct TetraP2::SubElem<4>::SNode<3> { enum { ID = 8 }; };

// Nodes of subelement 5
template <> template <> struct TetraP2::SubElem<5>::SNode<0> { enum { ID = 9 }; };
template <> template <> struct TetraP2::SubElem<5>::SNode<1> { enum { ID = 7 }; };
template <> template <> struct TetraP2::SubElem<5>::SNode<2> { enum { ID = 8 }; };
template <> template <> struct TetraP2::SubElem<5>::SNode<3> { enum { ID = 3 }; };

// Nodes of subelement 6
template <> template <> struct TetraP2::SubElem<6>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TetraP2::SubElem<6>::SNode<1> { enum { ID = 7 }; };
template <> template <> struct TetraP2::SubElem<6>::SNode<2> { enum { ID = 8 }; };
template <> template <> struct TetraP2::SubElem<6>::SNode<3> { enum { ID = 9 }; };

///////////////////////////////////////////////////////////////////////////////

/// shape function node 0
template <> template < typename COORD >
double TetraP2::ShapeF<0>::value () { return - lambda<COORD>() * ( 1. - 2.* lambda<COORD>()); }

/// shape function node 1
template <> template < typename COORD >
double TetraP2::ShapeF<1>::value () { return -COORD::xi() * (1. - 2.*COORD::xi()); }

/// shape function node 2
template <> template < typename COORD >
double TetraP2::ShapeF<2>::value () { return -COORD::eta() * (1. - 2.*COORD::eta()); }

/// shape function node 3
template <> template < typename COORD >
double TetraP2::ShapeF<3>::value () { return -COORD::zeta() * (1. - 2.*COORD::zeta()); }

/// shape function node 4
template <> template < typename COORD >
double TetraP2::ShapeF<4>::value () { return 4.*COORD::xi()*lambda<COORD>(); }

/// shape function node 5
template <> template < typename COORD >
double TetraP2::ShapeF<5>::value () { return 4.*COORD::xi()*COORD::eta(); }

/// shape function node 6
template <> template < typename COORD >
double TetraP2::ShapeF<6>::value () { return  4.*COORD::eta()*lambda<COORD>(); }

/// shape function node 7
template <> template < typename COORD >
double TetraP2::ShapeF<7>::value () { return 4.*COORD::xi()*COORD::zeta(); }

/// shape function node 8
template <> template < typename COORD >
double TetraP2::ShapeF<8>::value () { return 4.*COORD::eta()*COORD::zeta(); }

/// shape function node 9
template <> template < typename COORD >
double TetraP2::ShapeF<9>::value () { return 4.*COORD::zeta()*lambda<COORD>(); }

///////////////////////////////////////////////////////////////////////////////

template < unsigned int NELEM >
template < typename FSDATA >
void TetraP2::SubElem<NELEM>::compute_nodal_normals ( FSDATA& data )
{
  std::vector<RealVector>& nodal_normals = data.sub_elem[NELEM].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[NELEM].nodal_areas;

  typedef typename TetraP2::SubElem<NELEM>::template SNode<0> SNode0;
  typedef typename TetraP2::SubElem<NELEM>::template SNode<1> SNode1;
  typedef typename TetraP2::SubElem<NELEM>::template SNode<2> SNode2;
  typedef typename TetraP2::SubElem<NELEM>::template SNode<3> SNode3;

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

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_TetraP2_hh
