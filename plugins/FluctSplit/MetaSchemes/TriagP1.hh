#ifndef COOLFluiD_FluctSplit_TriagP1_hh
#define COOLFluiD_FluctSplit_TriagP1_hh

#include "FluctSplit/MetaSchemes/BaseSubElem.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

struct TriagP1
{
  enum { SHAPE      = CFGeoShape::TRIAG };
  enum { DIM        = DIM_2D };
  enum { ORDER      = CFPolyOrder::ORDER1 };
  enum { NBNODES    = 3 }; // (n^2+n)/2
  enum { NBSUBELEM  = Common::Power<ORDER,DIM>::value };
  enum { NBSUBNODES = 3 };

  static std::string getClassName() { return "TriagP1"; };
  static std::string getShapeName() { return CFGeoShape::Convert::to_str(CFGeoShape::TRIAG); };

  /// Class for each shape function.
  /// Template must be specialized.
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
template <> template <> struct TriagP1::SubElem<0>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TriagP1::SubElem<0>::SNode<1> { enum { ID = 1 }; };
template <> template <> struct TriagP1::SubElem<0>::SNode<2> { enum { ID = 2 }; };

///////////////////////////////////////////////////////////////////////////////

/// shape function node 0
/// \f[ N = 1 - \ksi - \eta - \zeta \f]
template <> template < typename COORD >
double TriagP1::ShapeF<0>::value ()
{ return 1.0 - COORD::xi() - COORD::eta(); }

/// shape function node 1
/// \f[ N = \eta \f]
template <> template < typename COORD >
double TriagP1::ShapeF<1>::value ()
{ return COORD::xi(); }

/// shape function node 2
/// \f[ N = \eta \f]
template <> template < typename COORD >
double TriagP1::ShapeF<2>::value ()
{ return COORD::eta(); }

///////////////////////////////////////////////////////////////////////////////

template <>
template < typename FSDATA >
void TriagP1::SubElem<0>::compute_nodal_normals ( FSDATA& data )
{
  std::vector<RealVector>& nodal_normals = data.sub_elem[0].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[0].nodal_areas;

  // nodes are switched to point normal inside
  compute_normal_vector( *data.nodes[2], *data.nodes[1], nodal_normals[0]);
  compute_normal_vector( *data.nodes[0], *data.nodes[2], nodal_normals[1]);
  compute_normal_vector( *data.nodes[1], *data.nodes[0], nodal_normals[2]);

  // compute unit normals and store nodal areas
  nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
  nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
  nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
}

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_TriagP1_hh
