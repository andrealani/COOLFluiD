#ifndef COOLFluiD_FluctSplit_TriagP2_hh
#define COOLFluiD_FluctSplit_TriagP2_hh

#include "FluctSplit/MetaSchemes/BaseSubElem.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluctSplit {

///////////////////////////////////////////////////////////////////////////////

struct TriagP2
{
  enum { SHAPE      = CFGeoShape::TRIAG };
  enum { DIM        = DIM_2D };
  enum { ORDER      = CFPolyOrder::ORDER2 };
  enum { NBNODES    = 6 }; // (n^2+n)/2
  enum { NBSUBELEM  = Common::Power<ORDER,DIM>::value };
  enum { NBSUBNODES = 3 };

  static std::string getClassName() { return "TriagP2"; };
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
template <> template <> struct TriagP2::SubElem<0>::SNode<0> { enum { ID = 0 }; };
template <> template <> struct TriagP2::SubElem<0>::SNode<1> { enum { ID = 3 }; };
template <> template <> struct TriagP2::SubElem<0>::SNode<2> { enum { ID = 5 }; };

// Nodes of subelement 1
template <> template <> struct TriagP2::SubElem<1>::SNode<0> { enum { ID = 3 }; };
template <> template <> struct TriagP2::SubElem<1>::SNode<1> { enum { ID = 1 }; };
template <> template <> struct TriagP2::SubElem<1>::SNode<2> { enum { ID = 4 }; };

// Nodes of subelement 2
template <> template <> struct TriagP2::SubElem<2>::SNode<0> { enum { ID = 5 }; };
template <> template <> struct TriagP2::SubElem<2>::SNode<1> { enum { ID = 4 }; };
template <> template <> struct TriagP2::SubElem<2>::SNode<2> { enum { ID = 2 }; };

// Nodes of subelement 3
template <> template <> struct TriagP2::SubElem<3>::SNode<0> { enum { ID = 5 }; };
template <> template <> struct TriagP2::SubElem<3>::SNode<1> { enum { ID = 3 }; };
template <> template <> struct TriagP2::SubElem<3>::SNode<2> { enum { ID = 4 }; };

///////////////////////////////////////////////////////////////////////////////

/// shape function node 0
/// \f[ N = \f]
template <> template < typename COORD >
double TriagP2::ShapeF<0>::value ()
{ return (1.0 - COORD::xi() - COORD::eta())*(2.*(1.0 - COORD::xi() - COORD::eta()) - 1.); }

/// shape function node 1
/// \f[ N =  \f]
template <> template < typename COORD >
double TriagP2::ShapeF<1>::value ()
{ return COORD::xi() * (2. * COORD::xi() - 1.); }

/// shape function node 2
/// \f[ N =  \f]
template <> template < typename COORD >
double TriagP2::ShapeF<2>::value ()
{ return COORD::eta() * (2. * COORD::eta() - 1.); }

/// shape function node 3
/// \f[ N =  \f]
template <> template < typename COORD >
double TriagP2::ShapeF<3>::value ()
{ return 4.*COORD::xi()*(1.0 - COORD::xi() - COORD::eta()); }

/// shape function node 4
/// \f[ N = \f]
template <> template < typename COORD >
double TriagP2::ShapeF<4>::value ()
{ return 4.*COORD::xi()*COORD::eta(); }

/// shape function node 5
/// \f[ N = \f]
template <> template < typename COORD >
double TriagP2::ShapeF<5>::value ()
{ return 4.*COORD::eta()*(1.0 - COORD::xi() - COORD::eta()); }

///////////////////////////////////////////////////////////////////////////////

// normals for subelem 0
template <>
template < typename FSDATA >
void TriagP2::SubElem<0>::compute_nodal_normals ( FSDATA& data )
{
  typedef TriagP2::SubElem<0> SElem;
  std::vector<RealVector>& nodal_normals = data.sub_elem[0].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[0].nodal_areas;
// CF_DEBUG_POINT;
//       CFout << "node [" << SElem::SNode<0>::ID << "]\n"; CFout << *data.nodes[SElem::SNode<0>::ID] << "\n";
//       CFout << "node [" << SElem::SNode<1>::ID << "]\n"; CFout << *data.nodes[SElem::SNode<1>::ID] << "\n";
//       CFout << "node [" << SElem::SNode<2>::ID << "]\n"; CFout << *data.nodes[SElem::SNode<2>::ID] << "\n";

      compute_normal_vector( *data.nodes[SElem::SNode<2>::ID], *data.nodes[SElem::SNode<1>::ID], nodal_normals[0]);
      compute_normal_vector( *data.nodes[SElem::SNode<0>::ID], *data.nodes[SElem::SNode<2>::ID], nodal_normals[1]);
      compute_normal_vector( *data.nodes[SElem::SNode<1>::ID], *data.nodes[SElem::SNode<0>::ID], nodal_normals[2]);

      // compute unit normals and store nodal areas
      nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
      nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
      nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
}

///////////////////////////////////////////////////////////////////////////////

// normals for subelem 1

template <>template < typename FSDATA >
void TriagP2::SubElem<1>::compute_nodal_normals ( FSDATA& data )
{
  typedef TriagP2::SubElem<1> SElem;
  std::vector<RealVector>& nodal_normals = data.sub_elem[1].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[1].nodal_areas;

// CF_DEBUG_POINT;
//       CFout << "node [" << SElem::SNode<0>::ID << "]\n"; CFout << *data.nodes[SElem::SNode<0>::ID] << "\n";
//       CFout << "node [" << SElem::SNode<1>::ID << "]\n"; CFout << *data.nodes[SElem::SNode<1>::ID] << "\n";
//       CFout << "node [" << SElem::SNode<2>::ID << "]\n"; CFout << *data.nodes[SElem::SNode<2>::ID] << "\n";

      compute_normal_vector( *data.nodes[SElem::SNode<2>::ID], *data.nodes[SElem::SNode<1>::ID], nodal_normals[0]);
      compute_normal_vector( *data.nodes[SElem::SNode<0>::ID], *data.nodes[SElem::SNode<2>::ID], nodal_normals[1]);
      compute_normal_vector( *data.nodes[SElem::SNode<1>::ID], *data.nodes[SElem::SNode<0>::ID], nodal_normals[2]);

      // compute unit normals and store nodal areas
      nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
      nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
      nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
}

///////////////////////////////////////////////////////////////////////////////

// normals for subelem 2
template <>
template < typename FSDATA >
void TriagP2::SubElem<2>::compute_nodal_normals ( FSDATA& data )
{
  typedef TriagP2::SubElem<2> SElem;
  std::vector<RealVector>& nodal_normals = data.sub_elem[2].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[2].nodal_areas;

      compute_normal_vector( *data.nodes[SElem::SNode<2>::ID], *data.nodes[SElem::SNode<1>::ID], nodal_normals[0]);
      compute_normal_vector( *data.nodes[SElem::SNode<0>::ID], *data.nodes[SElem::SNode<2>::ID], nodal_normals[1]);
      compute_normal_vector( *data.nodes[SElem::SNode<1>::ID], *data.nodes[SElem::SNode<0>::ID], nodal_normals[2]);

      // compute unit normals and store nodal areas
      nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
      nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
      nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
}

///////////////////////////////////////////////////////////////////////////////

// normals for subelem 3
template <>
template < typename FSDATA >
void TriagP2::SubElem<3>::compute_nodal_normals ( FSDATA& data )
{
  typedef TriagP2::SubElem<3> SElem;
  std::vector<RealVector>& nodal_normals = data.sub_elem[3].nodal_normals;
  std::vector<CFreal>&     nodal_areas   = data.sub_elem[3].nodal_areas;

      compute_normal_vector( *data.nodes[SElem::SNode<2>::ID], *data.nodes[SElem::SNode<1>::ID], nodal_normals[0]);
      compute_normal_vector( *data.nodes[SElem::SNode<0>::ID], *data.nodes[SElem::SNode<2>::ID], nodal_normals[1]);
      compute_normal_vector( *data.nodes[SElem::SNode<1>::ID], *data.nodes[SElem::SNode<0>::ID], nodal_normals[2]);

      // compute unit normals and store nodal areas
      nodal_areas[0] = nodal_normals[0].norm2(); nodal_normals[0] /= nodal_areas[0];
      nodal_areas[1] = nodal_normals[1].norm2(); nodal_normals[1] /= nodal_areas[1];
      nodal_areas[2] = nodal_normals[2].norm2(); nodal_normals[2] /= nodal_areas[2];
}

///////////////////////////////////////////////////////////////////////////////

  } // FluctSplit
} // COOLFluiD

#endif // COOLFluiD_FluctSplit_TriagP2_hh
