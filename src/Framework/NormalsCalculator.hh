// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NormalsCalculator_hh
#define COOLFluiD_Framework_NormalsCalculator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  class Node;
  class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a calculator for face normals of elements
/// @author Tiago Quintino
class Framework_API NormalsCalculator {
public:

  /// Default constructor without arguments.
  NormalsCalculator();

  /// Default destructor.
  ~NormalsCalculator();

  /// Compute the face normals for a Triangle
  /// with coordinates given in a matrix nodes x dimension.
  /// @param coord coordinates of the nodes
  /// @param normals the returned normals are passed in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (3x2)
  /// @pre normals must be sized (3x2)
  void computeTriagNormals(const RealMatrix& coord,
                                 RealMatrix& normals);

  /// Compute the unit face normals for a Triangle
  /// @param cell the GeometricEntity to calculate the normals
  /// @pre cell must be a Triag with P1 geometry
  void computeTriagUnitNormals(GeometricEntity& cell,
                               std::vector<RealVector>& unitNormals);

  /// Compute the face normals for a Quadrilateral
  /// with coordinates given in a matrix nodes x dimension.
  /// @param coord coordinates of the nodes
  /// @param normals the returned normals are passed in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (4x2)
  /// @pre normals must be sized (4x2)
  void computeQuadNormals(const RealMatrix& coord,
                                RealMatrix& normals);

  /// Compute the unit face normals for a Triangle
  /// @param cell the GeometricEntity to calculate the normals
  /// @pre cell must be a Quad with P1 geometry
  void computeQuadUnitNormals(GeometricEntity& cell,
                               std::vector<RealVector>& unitNormals);

  /// Compute the face normals for a tetrahedra.
  /// @param cell the GeometricEntity to calculate the normals
  /// @pre cell must be a Tetra with P1 geometry
  void computeTetraUnitNormals(GeometricEntity& cell,
                               std::vector<RealVector>& unitNormals);

  /// Compute the face normals for a tetrahedra
  /// with coordinates given in a matrix nodes x dimension.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates are passes in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (4x3)
  /// @pre normals must be sized (4x3)
  void computeTetraNormals(const RealMatrix& coord,
                                 RealMatrix& normals);

  /// Compute the face normals for a Pyramid.
  /// @param cell the GeometricEntity to calculate the normals
  /// @pre cell must be a Pyramid with P1 geometry
  void computePyramUnitNormals(GeometricEntity& cell,
                               std::vector<RealVector>& unitNormals);

  /// Compute the face normals for a pyramid with coordinates 
  /// given in a matrix nodes x dimension and check if they are outward.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates are passes in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (5x3)
  /// @pre normals must be sized (5x3)
  void computePyramOutwardNormals(const RealMatrix& coord,
				  RealMatrix& normals);
  
  /// Compute the face normals for a pyramid
  /// with coordinates given in a matrix nodes x dimension.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates are passes in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (5x3)
  /// @pre normals must be sized (5x3)
  void computePyramNormals(const RealMatrix& coord,
			   RealMatrix& normals);

  /// Compute the face normals for a Prism.
  /// @param cell the GeometricEntity to calculate the normals
  /// @pre cell must be a Prism with P1 geometry
  void computePrismUnitNormals(GeometricEntity& cell,
                               std::vector<RealVector>& unitNormals);

  /// Compute the face normals for a Prism
  /// with coordinates given in a matrix nodes x dimension.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates are passes in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (5x3)
  /// @pre normals must be sized (5x3)
  void computePrismNormals(const RealMatrix& coord,
                                 RealMatrix& normals);

  /// Compute the face normals for a Hexahedra
  /// @param cell the GeometricEntity to calculate the normals
  /// @pre cell must be a Hexa with P1 geometry
  void computeHexaUnitNormals(GeometricEntity& cell,
                               std::vector<RealVector>& unitNormals);

  /// Compute the face normals for a hexahedron
  /// with coordinates given in a matrix nodes x dimension.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates are passes in order
  ///                to avoid construction of another matrix
  /// @pre coord must be sized (8x3)
  /// @pre normals must be sized (6x3)
  void computeHexaNormals(const RealMatrix& coord,
                                RealMatrix& normals);


  /// Computes the normal to the triangle.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates
  void computeTriagNormal(RealVector& normal,
                          const RealMatrix& coord)
  {
    _vec1[XX] = coord(1,XX) - coord(0,XX);
    _vec1[YY] = coord(1,YY) - coord(0,YY);
    _vec1[ZZ] = coord(1,ZZ) - coord(0,ZZ);

    _vec2[XX] = coord(2,XX) - coord(1,XX);
    _vec2[YY] = coord(2,YY) - coord(1,YY);
    _vec2[ZZ] = coord(2,ZZ) - coord(1,ZZ);

    MathTools::MathFunctions::crossProd(_vec1,_vec2,normal);

    normal *= 0.5;
  }

  /// Computes the normal to a quad.
  /// @param coord coordinates of the nodes
  /// @param normals the returned coordinates
  void computeQuadNormal(RealVector& normal,
                         const RealMatrix& coord)
  {
    _vec1[XX] = coord(2,XX) - coord(0,XX);
    _vec1[YY] = coord(2,YY) - coord(0,YY);
    _vec1[ZZ] = coord(2,ZZ) - coord(0,ZZ);

    _vec2[XX] = coord(3,XX) - coord(1,XX);
    _vec2[YY] = coord(3,YY) - coord(1,YY);
    _vec2[ZZ] = coord(3,ZZ) - coord(1,ZZ);

    MathTools::MathFunctions::crossProd(_vec1,_vec2,normal);

    normal *= 0.5;
  }

protected:

  /// Computes the exterior normal to a line.
  /// Note that the order of indexes is important
  void computeFaceLineNormal(RealMatrix& normals,
                          const CFuint& n,
                          const RealMatrix& coord,
                          const CFuint& i,
                          const CFuint& j)
  {
    normals(n,XX) = coord(j,YY) - coord(i,YY);
    normals(n,YY) = coord(i,XX) - coord(j,XX);
  }

  /// Computes the exterior normal to a triangle.
  /// Note that the order of indexes is important
  /// Cycliv permutations are allowed:
  /// Fijk = Fjki = Fkij
  void computeFaceTriagNormal(RealMatrix& normals,
                          const CFuint& n,
                          const RealMatrix& coord,
                          const CFuint& i,
                          const CFuint& j,
                          const CFuint& k)
  {
    _vec1[XX] = coord(j,XX) - coord(i,XX);
    _vec1[YY] = coord(j,YY) - coord(i,YY);
    _vec1[ZZ] = coord(j,ZZ) - coord(i,ZZ);

    _vec2[XX] = coord(k,XX) - coord(j,XX);
    _vec2[YY] = coord(k,YY) - coord(j,YY);
    _vec2[ZZ] = coord(k,ZZ) - coord(j,ZZ);

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;

    normals(n,XX) = _vec3[XX];
    normals(n,YY) = _vec3[YY];
    normals(n,ZZ) = _vec3[ZZ];
  }

  /// Computes the exterior normal to a quadrilateral face.
  /// Note that the order of indexes is important
  /// Cyclic permutations are allowed:
  /// Fijkl = Fjkli = Fklij = Flijk
  void computeFaceQuadNormal(RealMatrix& normals,
                         const CFuint& n,
                         const RealMatrix& coord,
                         const CFuint& i,
                         const CFuint& j,
                         const CFuint& k,
                         const CFuint& l)
  {
    _vec1[XX] = coord(k,XX) - coord(i,XX);
    _vec1[YY] = coord(k,YY) - coord(i,YY);
    _vec1[ZZ] = coord(k,ZZ) - coord(i,ZZ);

    _vec2[XX] = coord(l,XX) - coord(j,XX);
    _vec2[YY] = coord(l,YY) - coord(j,YY);
    _vec2[ZZ] = coord(l,ZZ) - coord(j,ZZ);

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;

    normals(n,XX) = _vec3[XX];
    normals(n,YY) = _vec3[YY];
    normals(n,ZZ) = _vec3[ZZ];
  }


  /// Computes the exterior normal to a triangle.
  /// Note that the order of indexes is important
  /// Cycliv permutations are allowed:
  /// Fijk = Fjki = Fkij
  void computeFaceTriagOutwardNormal(const RealVector& centroid,
				     RealMatrix& normals,
				     const CFuint& n,
				     const RealMatrix& coord,
				     const CFuint& i,
				     const CFuint& j,
				     const CFuint& k)
  {
    _vec1[XX] = coord(j,XX) - coord(i,XX);
    _vec1[YY] = coord(j,YY) - coord(i,YY);
    _vec1[ZZ] = coord(j,ZZ) - coord(i,ZZ);

    _vec2[XX] = coord(k,XX) - coord(j,XX);
    _vec2[YY] = coord(k,YY) - coord(j,YY);
    _vec2[ZZ] = coord(k,ZZ) - coord(j,ZZ);

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;
    
    // check
    static RealVector cFace(3);
    static CFreal ov3 = 1./3.;
    for (CFuint dim = 0; dim < DIM_3D; ++dim) {
      cFace[dim] = ov3*(coord(i,dim) + coord(j,dim) + coord(k,dim));
    }
    _vec1 =  _vec3 - cFace;
    _vec2 =  centroid - cFace;
    if (MathTools::MathFunctions::innerProd(_vec1,_vec2) > 0.) {
      _vec3 *= -1.0;
    }
    
    normals(n,XX) = _vec3[XX];
    normals(n,YY) = _vec3[YY];
    normals(n,ZZ) = _vec3[ZZ];
  }

  /// Computes the exterior normal to a quadrilateral face.
  /// Note that the order of indexes is important
  /// Cyclic permutations are allowed:
  /// Fijkl = Fjkli = Fklij = Flijk
  void computeFaceQuadOutwardNormal(const RealVector& centroid,
				    RealMatrix& normals,
				    const CFuint& n,
				    const RealMatrix& coord,
				    const CFuint& i,
				    const CFuint& j,
				    const CFuint& k,
				    const CFuint& l)
  {
    _vec1[XX] = coord(k,XX) - coord(i,XX);
    _vec1[YY] = coord(k,YY) - coord(i,YY);
    _vec1[ZZ] = coord(k,ZZ) - coord(i,ZZ);

    _vec2[XX] = coord(l,XX) - coord(j,XX);
    _vec2[YY] = coord(l,YY) - coord(j,YY);
    _vec2[ZZ] = coord(l,ZZ) - coord(j,ZZ);

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;
    
    // check
    static RealVector cFace(3);
    for (CFuint dim = 0; dim < DIM_3D; ++dim) {
      cFace[dim] = 0.25*(coord(i,dim) + coord(j,dim) + coord(k,dim) + coord(l,dim));
    }
    
    _vec1 =  _vec3 - cFace;
    _vec2 =  centroid - cFace;
    if (MathTools::MathFunctions::innerProd(_vec1,_vec2) > 0.) {
      _vec3 *= -1.0;
    }
    
    normals(n,XX) = _vec3[XX];
    normals(n,YY) = _vec3[YY];
    normals(n,ZZ) = _vec3[ZZ];
  }
  
  /// Puts the nodes into a matrix of their coordinates to be used
  /// by the other functions
  void putNodesIntoCoord(const std::vector<Node*>* const nodes, RealMatrix& coordinates);

  /// Adimensionalize the 2d normals passed in by scaling them to unity lenght.
  void adimensionalize2d(const RealMatrix faceNormals, std::vector<RealVector>& unitNormals);

  /// Adimensionalize the 3d normals passed in by scaling them to unity lenght.
  void adimensionalize3d(const RealMatrix faceNormals, std::vector<RealVector>& unitNormals);

private: // data

  /// temporary matrix to compute the Quad normals
  RealMatrix _tmpNorm;

  /// temporary vector for computing Triag Face normals
  RealVector _vec1;

  /// temporary vector for computing Triag Face normals
  RealVector _vec2;

  /// temporary vector for computing Triag Face normals
  RealVector _vec3;

}; // end of class NormalsCalculator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_NormalsCalculator_hh
