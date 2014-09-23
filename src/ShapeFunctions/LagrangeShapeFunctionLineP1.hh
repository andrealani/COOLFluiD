// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionLineP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionLineP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/FaceJacobiansDeterminant.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

    class Node;

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a P1 (linear)
/// 1D element.
/// @author Andrea Lani
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionLineP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_1D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 2;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 2;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::LINE;
  }

  /// Gets the type of Interpolator
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::LAGRANGE;
  }

  /// Gets the Interpolator order
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER1;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {

    cf_assert (mappedCoords.size() == 2);

    mappedCoords[0][KSI] = 0. ;

    mappedCoords[1][KSI] = 1. ;


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
        const RealVector& mappedCoord, RealVector& shapeFunc)
  {
    cf_assert(mappedCoord.size() == 1);
    shapeFunc[0] = 0.5*(1.0 - mappedCoord[0]);
    shapeFunc[1] = 0.5*(1.0 + mappedCoord[0]);
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeGradientStates()");
  }

  /// Computes the normal to a face at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the dimensionality of the Face + 1)
  static void computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
                                                      const std::vector<Framework::Node*>& nodes,
                                                      std::vector<RealVector>& normal);

  /// Computes the normal to a given mapped coordinate plane, at the given mapped coordinates
  static void computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                            const std::vector<RealVector>& mappedCoord,
                                            const std::vector<Framework::Node*>& nodes,
                                            std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeMappedCoordPlaneNormal()");
  }

  /// Compute the Jacobian
  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobian()");
  }

  /// Compute the Jacobian
  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    const CFreal x1x0 = 0.5 * ((*nodes[1])[XX] - (*nodes[0])[XX]);
    const CFreal y1y0 = 0.5 * ((*nodes[1])[YY] - (*nodes[0])[YY]);

    const CFreal length = sqrt(x1x0*x1x0 + y1y0*y1y0);

    _tmpMat2D(XX,KSI) = x1x0;
    _tmpMat2D(XX,ETA) = y1y0;
    _tmpMat2D(YY,KSI) = -y1y0 / length;
    _tmpMat2D(YY,ETA) =  x1x0 / length;

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      jacob[ip] = _tmpMat2D;
    }
  }

  /// Compute the Jacobian
  static void computeJacobianPlus2D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    const CFreal x1x0 = 0.5 * ((*nodes[1])[XX] - (*nodes[0])[XX]);
    const CFreal y1y0 = 0.5 * ((*nodes[1])[YY] - (*nodes[0])[YY]);
    const CFreal z1z0 = 0.5 * ((*nodes[1])[ZZ] - (*nodes[0])[ZZ]);

    const CFreal length = sqrt(x1x0*x1x0 + y1y0*y1y0 + z1z0*z1z0);

    _tmpMat3D(XX,KSI) = x1x0;
    _tmpMat3D(XX,ETA) = 1.0 / length;
    _tmpMat3D(XX,ZTA) = y1y0;
    _tmpMat3D(YY,KSI) = y1y0;
    _tmpMat3D(YY,ETA) = 1.0 / length;
    _tmpMat3D(YY,ZTA) = z1z0;
    _tmpMat3D(ZZ,KSI) = z1z0;
    _tmpMat3D(ZZ,ETA) = 1.0 / length;
    _tmpMat3D(ZZ,ZTA) = x1x0;

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      jacob[ip] = _tmpMat3D;
    }
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminant(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    cf_assert(nodes.size() == getNbNodes());

    const CFreal x0 = (*nodes[0])[XX];
    const CFreal x1 = (*nodes[1])[XX];

    const CFreal jacob = 0.5 * (x1 - x0);

    detJacobian = jacob;
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates which have one more dimension than the
  /// GeometricEntity
  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    cf_assert(nodes.size() == getNbNodes());

    const CFreal x0 = (*nodes[0])[XX];
    const CFreal y0 = (*nodes[0])[YY];

    const CFreal x1 = (*nodes[1])[XX];
    const CFreal y1 = (*nodes[1])[YY];

    const CFreal dx = 0.5 * (x1 - x0);
    const CFreal dy = 0.5 * (y1 - y0);

    detJacobian = Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(dx,dy);
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates  which have two more dimension than the
  /// GeometricEntity
  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    cf_assert(nodes.size() == getNbNodes());

    const CFreal x0 = (*nodes[0])[XX];
    const CFreal y0 = (*nodes[0])[YY];
    const CFreal z0 = (*nodes[0])[ZZ];

    const CFreal x1 = (*nodes[1])[XX];
    const CFreal y1 = (*nodes[1])[YY];
    const CFreal z1 = (*nodes[1])[ZZ];

    const CFreal dx = 0.5 * (x1 - x0);
    const CFreal dy = 0.5 * (y1 - y0);
    const CFreal dz = 0.5 * (z1 - z0);

    detJacobian = Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(dx,dy,dz);
  }

  /// Get the name of this shape function
  static const std::string getName()
  {
    static const std::string name = "LagrangeLineP1";
    return name;
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

    const CFreal dx = (*nodes[0])[0] - (*nodes[1])[0];
    if (nodes[0]->size() == DIM_2D) {
      const CFreal dy = (*nodes[0])[1] - (*nodes[1])[1];
      return sqrt(dx*dx + dy*dy);
    }
    return std::abs(dx);
  }

  /// Get the centroid
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    RealVector centroid(*(nodes[0]));

    centroid += *(nodes[1]);
    centroid *= 0.5;
    return centroid;
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
    RealVector mappedCoord(coord.size() - DIM_1D);

    _tmpVec1 = (*nodes[0]);
    _tmpVec2 = (*nodes[1]);

//  CFout << "x0: " << _tmpVec1 << "\n";
//  CFout << "x1: " << _tmpVec2 << "\n";
//  CFout << "coord: " << coord << "\n";

    _tmpVec3 = _tmpVec2 - _tmpVec1;
    _tmpVec4 = coord - _tmpVec1;

    const CFreal normV = _tmpVec3.norm2();
    const CFreal normW = _tmpVec4.norm2();
    mappedCoord[XX] = normW/normV;
    //here we need to go from 0/1 to -1/1
    mappedCoord[XX] = (mappedCoord[XX]*2.) - 1.;
//  CFout << "Mapped Coord: " << mappedCoord << "\n";
    return mappedCoord;
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
    throw Common::ShouldNotBeHereException (FromHere(),getName()+"::computeAvgFaceNormals()");
  }

  /// Compute Face Average Normals
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeFaceNormals(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName()+"::computeFaceNormals()");
  }

  /// Compute Cell Average Normals
  /// @param nodes contains the nodes
  /// @return RealVector containing the average cell normal
  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    cf_assert(nodes.size() == 2);

    RealVector normal(DIM_2D);

    _tmpVec1 = (*nodes[0]);
    _tmpVec2 = (*nodes[1]);

    _vec = _tmpVec2 - _tmpVec1;
    normal[XX] = _vec[YY];
    normal[YY] = -_vec[XX];

    return normal;
  }

  /// Compute Face Average Normals
  /// @param nodes contains the nodes
  /// @return RealVector containing the cell normal
  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    /// Linear element -> same as computeAvgFaceNormals for all mappedCoord
    cf_assert(nodes.size() == 2);

    RealVector normal(DIM_2D);

    _tmpVec1 = (*nodes[0]);
    _tmpVec2 = (*nodes[1]);

    _vec = _tmpVec2 - _tmpVec1;
    normal[XX] = _vec[YY];
    normal[YY] = -_vec[XX];

    return normal;
  }

  /// Check if a point (defined with mapped coordinates) is inside an element
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 1);
    if((mappedCoord[0] < -1.) || (mappedCoord[0] > 1.))
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  /// Check if a point is inside an element
  static bool isInElement(const std::vector<Framework::Node*>& nodes, const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::isInElement()");
  }

private:

  /// Default constructor without arguments
  LagrangeShapeFunctionLineP1()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_2D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionLineP1();

private:

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// Temporary Vector for the computation of normals
  static RealVector _vec;

  /// Temporary Vector
  static RealVector _tmpVec1;

  /// Temporary Vector
  static RealVector _tmpVec2;

  /// Temporary Vector
  static RealVector _tmpVec3;

  /// Temporary Vector
  static RealVector _tmpVec4;

  /// Temporary Matrix
  static RealMatrix _tmpMat2D;

  /// Temporary Matrix
  static RealMatrix _tmpMat3D;

  /// Temporary Vector for the computation of normals
  static RealVector _vec3D;

  /// Vector of normals
  static std::vector<RealVector> _normals;

}; // end of class LagrangeShapeFunctionLineP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionLineP1_hh
