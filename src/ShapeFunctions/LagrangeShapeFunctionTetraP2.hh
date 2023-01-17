// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTetraP2_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTetraP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MatrixInverterT.hh"

#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/VolumeCalculator.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function
/// of the second order for a tetrahedron.
/// @author Thomas Wuilbaut
class ShapeFunctions_API LagrangeShapeFunctionTetraP2 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_3D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 10;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 4;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::TETRA;
  }

  /// Gets the type of Interpolator
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::LAGRANGE;
  }

  /// Gets the Interpolator order
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER2;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {

    cf_assert (mappedCoords.size() == 10);

    mappedCoords[0][KSI] = 0. ;
    mappedCoords[0][ETA] = 0. ;
    mappedCoords[0][ZTA] = 0. ;

    mappedCoords[1][KSI] = 1. ;
    mappedCoords[1][ETA] = 0. ;
    mappedCoords[1][ZTA] = 0. ;

    mappedCoords[2][KSI] = 0. ;
    mappedCoords[2][ETA] = 1. ;
    mappedCoords[2][ZTA] = 0. ;

    mappedCoords[3][KSI] = 0. ;
    mappedCoords[3][ETA] = 0. ;
    mappedCoords[3][ZTA] = 1. ;

    mappedCoords[4][KSI] = 0.5 ;
    mappedCoords[4][ETA] = 0. ;
    mappedCoords[4][ZTA] = 0. ;

    mappedCoords[5][KSI] = 0.5 ;
    mappedCoords[5][ETA] = 0.5 ;
    mappedCoords[5][ZTA] = 0. ;

    mappedCoords[6][KSI] = 0. ;
    mappedCoords[6][ETA] = 0.5 ;
    mappedCoords[6][ZTA] = 0. ;

    mappedCoords[7][KSI] = 0.5 ;
    mappedCoords[7][ETA] = 0. ;
    mappedCoords[7][ZTA] = 0.5 ;

    mappedCoords[8][KSI] = 0. ;
    mappedCoords[8][ETA] = 0.5 ;
    mappedCoords[8][ZTA] = 0.5 ;

    mappedCoords[9][KSI] = 0. ;
    mappedCoords[9][ETA] = 0. ;
    mappedCoords[9][ZTA] = 0.5 ;


  }

  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunctions(
          const std::vector<RealVector>& mappedCoord, std::vector<RealVector>& shapeFunc)
  {
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      computeShapeFunction(mappedCoord[ip],shapeFunc[ip]);
    }
  }

  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunction(
         const RealVector& mappedCoord,
               RealVector& shapeFunc)
  {
    cf_assert(mappedCoord.size() == DIM_3D);

    const CFreal lambda = 1.0 - mappedCoord.sum();
    const CFreal xi = mappedCoord[KSI];
    const CFreal eta = mappedCoord[ETA];
    const CFreal zeta = mappedCoord[ZTA];

    shapeFunc[0] = -lambda * (1. - 2.*lambda);
    shapeFunc[1] = -xi * (1. - 2.*xi);
    shapeFunc[2] = -eta * (1. - 2.*eta);
    shapeFunc[3] = -zeta * (1. - 2.*zeta);
    shapeFunc[4] = 4.*xi*lambda;
    shapeFunc[5] = 4.*xi*eta;
    shapeFunc[6] = 4.*eta*lambda;
    shapeFunc[7] = 4.*xi*zeta;
    shapeFunc[8] = 4.*eta*zeta;
    shapeFunc[9] = 4.*zeta*lambda;

  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    MathTools::MatrixInverterT<3> inverter;

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      inverter.invert(jacob[ip], m_invJ);

      RealMatrix& lgrad = grad[ip];
// CFout << "mappedCoord[ip]: " << mappedCoord[ip] <<"\n";
      const CFreal lambda = 1.0 - mappedCoord[ip].sum();
      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal dN0dxi = 1. - 4. * lambda;
      const CFreal dN1dxi = -1. + 4. * xi;
      const CFreal dN2dxi = 0.;
      const CFreal dN3dxi = 0.;
      const CFreal dN4dxi = 4.*(lambda - xi);
      const CFreal dN5dxi = 4. * eta;
      const CFreal dN6dxi = -4.* eta;
      const CFreal dN7dxi = 4. * zta;
      const CFreal dN8dxi = 0.;
      const CFreal dN9dxi = -4. * zta;

      const CFreal dN0deta = 1. - 4. * lambda;
      const CFreal dN1deta = 0.;
      const CFreal dN2deta = -1. + 4. * eta;
      const CFreal dN3deta = 0.;
      const CFreal dN4deta = -4. * xi;
      const CFreal dN5deta = 4. * xi;
      const CFreal dN6deta = 4.*(lambda - eta);
      const CFreal dN7deta = 0.;
      const CFreal dN8deta = 4. * zta;
      const CFreal dN9deta = -4. * zta;

      const CFreal dN0dzta = 1. - 4. * lambda;
      const CFreal dN1dzta = 0.;
      const CFreal dN2dzta = 0.;
      const CFreal dN3dzta = -1. + 4. * zta;
      const CFreal dN4dzta = -4. * xi;
      const CFreal dN5dzta = 0.;
      const CFreal dN6dzta = -4. * eta;
      const CFreal dN7dzta = 4. * xi;
      const CFreal dN8dzta = 4. * eta;
      const CFreal dN9dzta = 4. * (lambda - zta);

      const CFreal JXX = m_invJ(0,0);
      const CFreal JXY = m_invJ(0,1);
      const CFreal JXZ = m_invJ(0,2);
      const CFreal JYX = m_invJ(1,0);
      const CFreal JYY = m_invJ(1,1);
      const CFreal JYZ = m_invJ(1,2);
      const CFreal JZX = m_invJ(2,0);
      const CFreal JZY = m_invJ(2,1);
      const CFreal JZZ = m_invJ(2,2);

      lgrad(0,XX) = dN0dxi * JXX + dN0deta * JXY + dN0dzta * JXZ;
      lgrad(0,YY) = dN0dxi * JYX + dN0deta * JYY + dN0dzta * JYZ;
      lgrad(0,ZZ) = dN0dxi * JZX + dN0deta * JZY + dN0dzta * JZZ;

/*CFout << "grad(0,XX): " << lgrad(0,XX) <<"\n";
CFout << "dN0dxi: " << dN0dxi <<"\n";
CFout << "dN0deta: " << dN0deta <<"\n";
CFout << "dN0dzta: " << dN0dzta <<"\n";
CFout << "JXX: " << JXX <<"\n";
CFout << "JXY: " << JXY <<"\n";
CFout << "JXZ: " << JXZ <<"\n";*/
      lgrad(1,XX) = dN1dxi * JXX + dN1deta * JXY + dN1dzta * JXZ;
      lgrad(1,YY) = dN1dxi * JYX + dN1deta * JYY + dN1dzta * JYZ;
      lgrad(1,ZZ) = dN1dxi * JZX + dN1deta * JZY + dN1dzta * JZZ;

      lgrad(2,XX) = dN2dxi * JXX + dN2deta * JXY + dN2dzta * JXZ;
      lgrad(2,YY) = dN2dxi * JYX + dN2deta * JYY + dN2dzta * JYZ;
      lgrad(2,ZZ) = dN2dxi * JZX + dN2deta * JZY + dN2dzta * JZZ;

      lgrad(3,XX) = dN3dxi * JXX + dN3deta * JXY + dN3dzta * JXZ;
      lgrad(3,YY) = dN3dxi * JYX + dN3deta * JYY + dN3dzta * JYZ;
      lgrad(3,ZZ) = dN3dxi * JZX + dN3deta * JZY + dN3dzta * JZZ;

      lgrad(4,XX) = dN4dxi * JXX + dN4deta * JXY + dN4dzta * JXZ;
      lgrad(4,YY) = dN4dxi * JYX + dN4deta * JYY + dN4dzta * JYZ;
      lgrad(4,ZZ) = dN4dxi * JZX + dN4deta * JZY + dN4dzta * JZZ;

      lgrad(5,XX) = dN5dxi * JXX + dN5deta * JXY + dN5dzta * JXZ;
      lgrad(5,YY) = dN5dxi * JYX + dN5deta * JYY + dN5dzta * JYZ;
      lgrad(5,ZZ) = dN5dxi * JZX + dN5deta * JZY + dN5dzta * JZZ;

      lgrad(6,XX) = dN6dxi * JXX + dN6deta * JXY + dN6dzta * JXZ;
      lgrad(6,YY) = dN6dxi * JYX + dN6deta * JYY + dN6dzta * JYZ;
      lgrad(6,ZZ) = dN6dxi * JZX + dN6deta * JZY + dN6dzta * JZZ;

      lgrad(7,XX) = dN7dxi * JXX + dN7deta * JXY + dN7dzta * JXZ;
      lgrad(7,YY) = dN7dxi * JYX + dN7deta * JYY + dN7dzta * JYZ;
      lgrad(7,ZZ) = dN7dxi * JZX + dN7deta * JZY + dN7dzta * JZZ;

      lgrad(8,XX) = dN8dxi * JXX + dN8deta * JXY + dN8dzta * JXZ;
      lgrad(8,YY) = dN8dxi * JYX + dN8deta * JYY + dN8dzta * JYZ;
      lgrad(8,ZZ) = dN8dxi * JZX + dN8deta * JZY + dN8dzta * JZZ;

      lgrad(9,XX) = dN9dxi * JXX + dN9deta * JXY + dN9dzta * JXZ;
      lgrad(9,YY) = dN9dxi * JYX + dN9deta * JYY + dN9dzta * JYZ;
      lgrad(9,ZZ) = dN9dxi * JZX + dN9deta * JZY + dN9dzta * JZZ;
    }
  }

  /// Computes the normal to a face at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the dimensionality of the Face + 1)
  static void computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
      const std::vector<Framework::Node*>& nodes,
      std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeFaceJacobDetVectorAtMappedCoords()");
  }

  /// Computes the normal to a given mapped coordinate plane, at the given mapped coordinates
  static void computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                            const std::vector<RealVector>& mappedCoord,
                                            const std::vector<Framework::Node*>& nodes,
                                            std::vector<RealVector>& normal);

  /// Compute the Jacobian
  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
                              std::vector<RealMatrix>& jacob);

  /// Compute the Jacobian
  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianPlus1D()");
  }

  /// Compute the Jacobian
  static void computeJacobianPlus2D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianPlus2D()");
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminant(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
                                         std::valarray<CFreal>& detJacobian)
    {
      /// @note implemented in a not very efficient way
      cf_assert(nodes.size() == getNbNodes());

      const CFuint nbrPnts = mappedCoord.size();
      cf_assert(detJacobian.size() == nbrPnts);

      std::vector< RealMatrix > pntJacob(nbrPnts,RealMatrix(3,3));
      
      computeJacobian(nodes,mappedCoord,pntJacob);
      
      for (CFuint ip = 0; ip < nbrPnts; ++ip)
      {
        detJacobian[ip] = pntJacob[ip].determ3();
      }
    }


  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobianDeterminantPlus1D()");
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianDeterminantPlus2D()");
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  /// @todo fix this function
  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian)
  {

    throw Common::NotImplementedException (FromHere(),getName()
          + "::computeFaceJacobianDeterminant()");
  }

  /// Get the name of this shape function
  static const std::string getName()
  {
    return "LagrangeTetraP2";
  }

  /// Get the ID for the solution integration
  static Framework::InterpolatorID getInterpolatorID()
  {
    return _interpolatorID;
  }

  /// Get the ID for the solution integration
  static void setInterpolatorID(const Framework::InterpolatorID& id)
  {
    _interpolatorID = id;
  }

  /// Get the volume
  static CFreal computeVolume(const std::vector<Framework::Node*>& nodes)
  {
    using namespace std;
    using namespace COOLFluiD::Framework;


    static RealMatrix coord(4,DIM_3D);
    static VolumeCalculator volumeCalc;

    for (CFuint iNode = 0; iNode < 4; ++iNode) {
      for (CFuint iCoord = 0; iCoord < DIM_3D; ++iCoord) {
        coord(iNode, iCoord) = (*nodes[iNode])[iCoord];
      }
    }
    return volumeCalc.calculateTetraVolume(coord);
  }

  /// Get the centroid
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeCentroid()");

  }

  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinates(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
      throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinates()");
  }


  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinatesPlus1D(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinatesPlus1D() ");
  }

  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinatesPlus2D(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinatesPlus2D()");
  }

  /// Compute Face Average Normals
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeAvgFaceNormals(const std::vector<Framework::Node*>& nodes)
  {
    cf_assert(nodes.size() == 4);

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);

    // Face 013
    computeFaceTriagNormal(1,nodes,0,1,3);

    // Face 123
    computeFaceTriagNormal(2,nodes,1,2,3);

    // Face 032
    computeFaceTriagNormal(3,nodes,0,3,2);

    return _normals;
  }

  /// Compute Face Normals
  /// @param mappedCoord contains the coordinates of the location at which the face normal should be computed
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    /// Linear element -> same as computeAvgFaceNormals for all mappedCoord
    cf_assert(nodes.size() == 4);

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);

    // Face 013
    computeFaceTriagNormal(1,nodes,0,1,3);

    // Face 123
    computeFaceTriagNormal(2,nodes,1,2,3);

    // Face 032
    computeFaceTriagNormal(3,nodes,0,3,2);

    return _normals;
  }

  /// Compute Cell Average Normal
  /// @param nodes contains the nodes
  /// @return RealVector containing the average cell normal
  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName()+"::computeAvgCellNormal()");
  }

  /// Compute Cell Normal at a given mapped coord
  /// @param nodes contains the nodes
  /// @return RealVector containing the cell normal
  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName()+"::computeCellNormal()");
  }

  /// Check if a point (defined with mapped coordinates) is inside an element
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 3);
    if( (mappedCoord[0] >= 0.) &&
        (mappedCoord[1] >= 0.) &&
        (mappedCoord[2] >= 0.) &&
        (mappedCoord.sum() <= 1.))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  /// Check if a point is inside an element
  static bool isInElement(const std::vector<Framework::Node*>& nodes, const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::isInElement()");
  }

private:

  /// Constructor without arguments
  LagrangeShapeFunctionTetraP2()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionTetraP2();

  /// Computes the normal to the triangle.
  /// Note that the order of indexes is important
  /// Cycliv permutations are allowed:
  /// Fijk = Fjki = Fkij
  static void computeFaceTriagNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                          const CFuint& i,
                          const CFuint& j,
                          const CFuint& k)
  {
    _vec1[XX] = (*nodes[j])[XX] - (*nodes[i])[XX];
    _vec1[YY] = (*nodes[j])[YY] - (*nodes[i])[YY];
    _vec1[ZZ] = (*nodes[j])[ZZ] - (*nodes[i])[ZZ];

    _vec2[XX] = (*nodes[k])[XX] - (*nodes[j])[XX];
    _vec2[YY] = (*nodes[k])[YY] - (*nodes[j])[YY];
    _vec2[ZZ] = (*nodes[k])[ZZ] - (*nodes[j])[ZZ];

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;

    _normals[n][XX] = _vec3[XX];
    _normals[n][YY] = _vec3[YY];
    _normals[n][ZZ] = _vec3[ZZ];
  }

private:

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// Temporary matrix for the holding the 2D jacobian inverse
  static RealMatrix m_invJ;

  /// Temporary Vector for the computation of normals
  static RealVector _vec1;
  static RealVector _vec2;
  static RealVector _vec3;
  static RealVector _vec4;
    
  /// gradients of shape functions
  static std::vector< RealVector > _gradShapFunc;

  /// Vector of normals
  static RealVector m_mappedCoord;

  /// Temporary Matrices
  static RealMatrix m_matrix1;
  static RealMatrix m_matrix2;

  /// Vector of normals
  static std::vector<RealVector> _normals;

}; // end of class LagrangeShapeFunctionTetraP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionTetraP2.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTetraP2_hh
