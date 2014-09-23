#ifndef COOLFluiD_Numerics_FluctSplit_P2Normal_hh
#define COOLFluiD_Numerics_FluctSplit_P2Normal_hh


#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class computes normals at arbitrary point of boundary of a P2P2 triangle
/// @author Martin Vymazal
class P2Normal {

public:

  /// Default constructor without arguments
  P2Normal();

  /// Default destructor
  ~P2Normal();

  /// Returns normal perpendicular to a triangle face at given point in reference space
  /// @param nodes nodes of P2P2 triangle
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

  CFuint m_workFaceIdx;
  RealVector m_weights;
  RealVector m_qpPos; //Position of quadrature points
                      //Integration domain <-0.25;0.25> considered


  //Shift coordinates to have interval
  // a) <0.0 ; 0.5>
  // b) <0.5 ; 1.0>
  //The reference triangle has edge defined for xi = <0.0; 1.0>
  //We need two normals for this edge for P2 triangle, therefore two
  //intervals

  RealVector m_shift;

}; // end of class P2Normal

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




#endif // COOLFluiD_Numerics_FluctSplit_P2Normal_hh
