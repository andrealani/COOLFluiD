#ifndef COOLFluiD_Numerics_FluctSplit_P3Normal_hh
#define COOLFluiD_Numerics_FluctSplit_P3Normal_hh


#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class computes normals at arbitrary point of boundary of a P3P3 (sub)triangle
/// @author Martin Vymazal
class P3Normal {

public:

  /// Default constructor without arguments
  P3Normal();

  /// Default destructor
  ~P3Normal();

  /// Returns normal perpendicular to a triangle face at given point in reference space
  /// @param nodes nodes of P3P3 triangle
  /// @param faceIdx defines at which of the nine subfaces the integration point lies
  /// @param xi,eta coordinates of the integration point in reference triangle
  /// @return coordinates of the normal
void ComputeNormal(const std::vector<Framework::Node*>& nodes, const CFuint faceIdx, const CFreal xi, const CFreal eta, RealVector& normal);


  /// Computes face lenghts of a curved P2P2 triangle
  /// @param nodes nodes of P2P2 triangle
//  void SetFaceLengths(const std::vector<Framework::Node*>& nodes);


  /// Returns normal perpendicular to boundary face at given point in reference space
  /// @param nodes nodes of P2P2 triangle
  /// @param faceIdx defines at which of the nine subfaces the integration point lies
  /// @param xi,eta coordinates of the integration point in reference triangle
  /// @return coordinates of the normal
void ComputeBNormal(const std::vector<Framework::Node*>& nodes, const CFreal xi, RealVector& normal);


private:

  static const CFreal m_XI[10];
  static const CFreal m_ETA[10];
  static const CFreal m_dLdt_base[3][3]; //derivatives of L0, L1, L2 along the subfaces 0,1,2 of a P3P3 triangle
                                         //first row contains derivatives of L0, L1, L2 along face 0
  static const CFuint m_nVert;
  static const CFuint m_nSubFaces;

  CFreal m_dndt[10];
  CFreal m_w[10];
  CFreal m_L[3];
  CFreal const *m_dLdt;

  CFreal m_dN0dxi, m_dN1dxi, m_dN2dxi, m_dN3dxi;


}; // end of class P3Normal

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




#endif // COOLFluiD_Numerics_FluctSplit_P3Normal_hh
