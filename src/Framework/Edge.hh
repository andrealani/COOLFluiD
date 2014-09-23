// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Edge_hh
#define COOLFluiD_Framework_Edge_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Edge. It is parameterized with the concrete
/// type of ShapeFunction, tecnique that allows to make inlinable
/// all the public methods of the shape function when used by this class.
/// The difference between a function call and a virtual function call
/// is smaller in terms of performance than the difference between an
/// inline method and the not inlined version of the same method.
/// The design choice of putting the ShapeFunction as a template parameter of
/// the Edge is not less flexible that the one of putting a pointer to the base
/// class of the shape functions, but is surely more efficient.
/// For more info about this tecnique, you can look at "Efficient C++" by
/// Bulka and Mayhew
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename GEO_SHAPE_FUNCTION,
          typename SOL_SHAPE_FUNCTION>
class Edge : public GeometricEntity {
public:

  /// Constructor
  /// The area the Edge is computed and stored
  /// @param states  list of the states in the Edge
  /// @param nodes list of the nodes in the Edge
  Edge(const std::vector<State*>& states,
      const std::vector<Node*>& nodes) :
    GeometricEntity(states, nodes)
  {
    cf_assert(GEO_SHAPE_FUNCTION::getShape() == SOL_SHAPE_FUNCTION::getShape());
  }

  /// Default destructor
  ~Edge() {}

  /// Computes the volume, area or length
  /// @pre Are considered volume of a geometric entity
  ///      in 3D: cell volume, Edge area, edge length;
  ///      in 2D: cell area, Edge length;
  ///      in 1D: cell length.
  CFreal computeVolume()
  {
    return GEO_SHAPE_FUNCTION::computeVolume(_nodes);
  }

  /// Computes the centroid
  RealVector computeCentroid()
  {
    return GEO_SHAPE_FUNCTION::computeCentroid(_nodes);
  }

  /// Get the inheritant dimensionality of the Face
  CFuint getDimensionality() const
  {
    cf_assert(GEO_SHAPE_FUNCTION::getDimensionality() == SOL_SHAPE_FUNCTION::getDimensionality());
    return GEO_SHAPE_FUNCTION::getDimensionality();
  }

  /// Gets the CFGeoEnt::Type of the Geometric Entity
  /// @post CFGeoEnt::Type is CFGeoEnt::Type::EDGE
  static CFGeoEnt::Type getGeomType()
  {
    return CFGeoEnt::EDGE;
  }

  /// Gets the CFGeoEnt::Type of the Geometric Entity
  /// @post CFGeoEnt::Type is CFGeoEnt::Type::EDGE
  virtual CFGeoEnt::Type getGeometryType() const
  {
    return Edge<GEO_SHAPE_FUNCTION,SOL_SHAPE_FUNCTION>::getGeomType();
  }

  /// Gets the type of CFGeoShape::Type
  const CFGeoShape::Type getShape() const
  {
    cf_assert(GEO_SHAPE_FUNCTION::getShape() == SOL_SHAPE_FUNCTION::getShape());
    return GEO_SHAPE_FUNCTION::getShape();
  }

  /// Gets the number of nodes of the geometric shape function
  CFuint getNbNodesGeometryShapeFunction() const
  {
    return GEO_SHAPE_FUNCTION::getNbNodes();
  }

  /// Gets the number of nodes of the solution shape function
  CFuint getNbNodesSolutionShapeFunction() const
  {
    return SOL_SHAPE_FUNCTION::getNbNodes();
  }

  /// Gets the ID for the solution interpolation
  InterpolatorID getSolInterpolatorID() const
  {
    return SOL_SHAPE_FUNCTION::getInterpolatorID();
  }

  /// Gets the ID for the geometric interpolation
  InterpolatorID getGeomInterpolatorID() const
  {
    return GEO_SHAPE_FUNCTION::getInterpolatorID();
  }

  /// Gets the type of geometry shape function
  const CFPolyForm::Type getGeometryShapeFunctionType() const
  {
    return GEO_SHAPE_FUNCTION::getInterpolatorType();
  }

  /// Gets the geometry shape function order
  const CFPolyOrder::Type getGeometryShapeFunctionOrder() const
  {
    return GEO_SHAPE_FUNCTION::getInterpolatorOrder();
  }

  /// Gets the type of solution shape function
  const CFPolyForm::Type getSolutionShapeFunctionType() const
  {
    return SOL_SHAPE_FUNCTION::getInterpolatorType();
  }

  /// Gets the solution shape function  order
  const CFPolyOrder::Type getSolutionShapeFunctionOrder() const
  {
    return GEO_SHAPE_FUNCTION::getInterpolatorOrder();
  }

  /// Computes the gradients of the shapefunctions at given coordinates
  std::vector<RealMatrix> computeGeoShapeFunctionGradients(const std::vector<RealVector>& coord)
  {
    throw Common::NotImplementedException
        (FromHere(), "Edge::computeGeoShapeFunctionGradients()");
  }

  /// Computes the gradients of the shapefunctions at given coordinates
  std::vector<RealMatrix> computeSolutionShapeFunctionGradients(const std::vector<RealVector>& coord)
  {
  // Compute Jacobian
  CFuint dim = PhysicalModelStack::getActive()->getDim();
  CFuint nbNodes = coord.size();
  std::vector<RealMatrix> jacob(nbNodes);
  std::vector<RealMatrix> grad(nbNodes);
  std::vector<RealVector> mappedCoord(nbNodes);
  for(CFuint i=0; i<nbNodes ; ++i){
    jacob[i].resize(dim,dim);
    grad[i].resize(nbNodes,dim);
    mappedCoord[i].resize(coord[0].size());
  }

  if(dim == SOL_SHAPE_FUNCTION::getDimensionality()) {
    for (CFuint i=0; i<coord.size();i++){
      mappedCoord[i].resize(dim);
      mappedCoord[i] = SOL_SHAPE_FUNCTION::computeMappedCoordinates(coord[i],_nodes);
      }
    SOL_SHAPE_FUNCTION::computeJacobian(_nodes, mappedCoord, jacob);
  }
  else {
    if(dim - DIM_1D == SOL_SHAPE_FUNCTION::getDimensionality()) {
      for (CFuint i=0; i<coord.size();i++){
        mappedCoord[i].resize(dim - DIM_1D);
        mappedCoord[i] = SOL_SHAPE_FUNCTION::computeMappedCoordinatesPlus1D(coord[i],_nodes);
      }
      SOL_SHAPE_FUNCTION::computeJacobianPlus1D(_nodes, mappedCoord, jacob);
    }
    if(dim - DIM_2D == SOL_SHAPE_FUNCTION::getDimensionality()) {
      for (CFuint i=0; i<coord.size();i++){
        mappedCoord[i] = SOL_SHAPE_FUNCTION::computeMappedCoordinatesPlus2D(coord[i],_nodes);
      }
      SOL_SHAPE_FUNCTION::computeJacobianPlus2D(_nodes, mappedCoord, jacob);
    }
  }

  // Compute the Gradients of the shape function at the given coord
  SOL_SHAPE_FUNCTION::computeGradientStates(jacob,mappedCoord,grad);

  return grad;
  }

  /// Computes the gradients of the solution shapefunctions at given coordinates
  std::vector<RealMatrix> computeSolutionShapeFunctionGradientsInMappedCoordinates(const std::vector<RealVector>& coord) const
  {
    throw Common::NotImplementedException
        (FromHere(), "Edge::computeSolutionShapeFunctionGradientsInMappedCoordinates()");
  }

  /// Computes the Jacobian matrix of the geometric shapefunctions at given coordinates
  std::vector<RealMatrix> computeGeometricShapeFunctionJacobianMatrix(const std::vector<RealVector>& coord) const
  {
    throw Common::NotImplementedException
        (FromHere(), "Edge::computeGeometricShapeFunctionJacobianMatrix()");
  }

  /// Computes the Jacobian determinant of the geometric shapefunctions at given mapped coordinates
  std::valarray<CFreal> computeGeometricShapeFunctionJacobianDeterminant(const std::vector<RealVector>& mappedCoord) const
  {
    throw Common::NotImplementedException
        (FromHere(), "Edge::computeGeometricShapeFunctionJacobianMatrix()");
  }

  /// Computes the mapped coordinates from the coordinates
  RealVector computeMappedCoordFromCoord(const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), "Edge::computeMappedCoordFromCoord()");
  }

  /// Computes the value of the shapefunctions at a given set of coordinates
  RealVector computeShapeFunctionAtCoord(const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), "Edge::computeShapeFunctionAtCoord()");
  }

  /// Computes the value of the shapefunctions at a given set of coordinates
  RealVector computeGeoShapeFunctionAtCoord(const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), "Edge::computeShapeFunctionAtCoord()");
  }

  /// Computes the coordinates from the mapped coordinates
  RealVector computeCoordFromMappedCoord(const RealVector& coord)
  {
    throw Common::NotImplementedException (FromHere(),"Edge::computeCoordFromMappedCoord()");
  }

  /// Computes the shape function from the mapped coordinates
  RealVector computeShapeFunctionAtMappedCoord(const RealVector& mappedCoord)
  {
    const CFuint nbStates = _states.size();
    RealVector shapeFunc(nbStates);

    SOL_SHAPE_FUNCTION::computeShapeFunction(mappedCoord,shapeFunc);

    return shapeFunc;
  }

  /// Computes the shape function from the mapped coordinates
  RealVector computeGeoShapeFunctionAtMappedCoord(const RealVector& mappedCoord)
  {
    const CFuint nbNodes = _nodes.size();
    RealVector shapeFunc(nbNodes);

    GEO_SHAPE_FUNCTION::computeShapeFunction(mappedCoord,shapeFunc);

    return shapeFunc;
  }

  /// Computes the average Face Normals
  /// (Normals have the dimensionality of the Cell)
  std::vector<RealVector> computeAvgFaceNormals()
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::computeAvgFaceNormals()");
  }

  /// Computes the Face Normals at the given coordinates
  /// (Normals have the dimensionality of the Cell)
  std::vector<RealVector> computeFaceNormals(const RealVector& coord)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::computeFaceNormals()");
  }

  /// Computes the normal to a given mapped coordinate plane, at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the same dimensionality as the Cell)
  std::vector< RealVector >
      computeMappedCoordPlaneNormalAtMappedCoords(const std::vector< CFuint >& planeIdx,
                                                  const std::vector< RealVector >& coord)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::computeMappedCoordPlaneNormalAtMappedCoords()");
  }

  /// Computes the normal to a face at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the dimensionality of the Face + 1)
  std::vector< RealVector >
      computeFaceJacobDetVectorAtMappedCoords(const std::vector< RealVector >& coord)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::computeFaceJacobDetVectorAtMappedCoords()");
  }

  /// Computes the average normals of the Cell
  /// (Normal has the dimensionality as the Cell + 1)
  RealVector computeAvgCellNormal()
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::computeAvgCellNormal()");
  }

  /// Computes the average normals of the Cell at given coordinates
  /// (Normal has the dimensionality as the Cell + 1)
  RealVector computeCellNormal(const RealVector& coord)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::computeCellNormal()");
  }

  /// Check if a point (defined with mapped coordinates) is inside an element
  bool isInMappedElement(const RealVector& mappedCoord)
  {
    return GEO_SHAPE_FUNCTION::isInMappedElement(mappedCoord);
  }

  /// Check if a point is inside an element
  bool isInElement(const RealVector& coord)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Edge::isInElement()");
  }

private: // data

}; // end of class Edge

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Edge_hh
