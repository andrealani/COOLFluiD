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

#include "NormalsCalculator.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/MathChecks.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

NormalsCalculator::NormalsCalculator() :
  _tmpNorm(1,DIM_3D),
  _vec1(DIM_3D),
  _vec2(DIM_3D),
  _vec3(DIM_3D)
{
}

//////////////////////////////////////////////////////////////////////////////

NormalsCalculator::~NormalsCalculator()
{
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeTetraNormals(const RealMatrix& coord,
                                                  RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 4);
  cf_assert(coord.nbCols() == DIM_3D);

  cf_assert(normals.nbRows() == 4);
  cf_assert(normals.nbCols() == DIM_3D);

  // Face 021
  computeFaceTriagNormal(normals,0,coord,0,2,1);

  // Face 013
  computeFaceTriagNormal(normals,1,coord,0,1,3);

  // Face 123
  computeFaceTriagNormal(normals,2,coord,1,2,3);

  // Face 032
  computeFaceTriagNormal(normals,3,coord,0,3,2);

}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computePyramOutwardNormals(const RealMatrix& coord,
						   RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 5);
  cf_assert(coord.nbCols() == DIM_3D);
  
  cf_assert(normals.nbRows() == 5);
  cf_assert(normals.nbCols() == DIM_3D);
  
  // centroid of the cell
  static RealVector cVolume(3);
  for (CFuint dim = 0; dim < 3; ++dim)  {
    cVolume[dim] = 0.2*(coord(0,dim) + coord(1,dim) +
			coord(2,dim) + coord(3,dim) +
			coord(4,dim));
  }
  
  // Face 0321
  computeFaceQuadOutwardNormal(cVolume, normals,0,coord,0,3,2,1);
  
  // Face 014
  computeFaceTriagOutwardNormal(cVolume,normals,1,coord,0,1,4);
  
  // Face 124
  computeFaceTriagOutwardNormal(cVolume,normals,2,coord,1,2,4);
  
  // Face 234
  computeFaceTriagOutwardNormal(cVolume,normals,3,coord,2,3,4);
  
  // Face 043
  computeFaceTriagOutwardNormal(cVolume,normals,4,coord,0,4,3);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computePyramNormals(const RealMatrix& coord,
					    RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 5);
  cf_assert(coord.nbCols() == DIM_3D);

  cf_assert(normals.nbRows() == 5);
  cf_assert(normals.nbCols() == DIM_3D);

  // Face 0321
  computeFaceQuadNormal(normals,0,coord,0,3,2,1);

  // Face 014
  computeFaceTriagNormal(normals,1,coord,0,1,4);

  // Face 124
  computeFaceTriagNormal(normals,2,coord,1,2,4);

  // Face 234
  computeFaceTriagNormal(normals,3,coord,2,3,4);

  // Face 043
  computeFaceTriagNormal(normals,4,coord,0,4,3);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computePrismNormals(const RealMatrix& coord,
                                                  RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 6);
  cf_assert(coord.nbCols() == DIM_3D);

  cf_assert(normals.nbRows() == 5);
  cf_assert(normals.nbCols() == DIM_3D);

  // Face 021
  computeFaceTriagNormal(normals,0,coord,0,2,1);

  // Face 345
  computeFaceTriagNormal(normals,1,coord,3,4,5);

  // Face 0143
  computeFaceQuadNormal(normals,2,coord,0,1,4,3);

  // Face 1254
  computeFaceQuadNormal(normals,3,coord,1,2,5,4);

  // Face 0352
  computeFaceQuadNormal(normals,4,coord,0,3,5,2);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeHexaNormals(const RealMatrix& coord,
                                                 RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 8);
  cf_assert(coord.nbCols() == DIM_3D);

  cf_assert(normals.nbRows() == 6);
  cf_assert(normals.nbCols() == 3);

  // Face 0321
  computeFaceQuadNormal(normals,0,coord,0,3,2,1);

  // Face 4567
  computeFaceQuadNormal(normals,1,coord,4,5,6,7);

  // Face 0154
  computeFaceQuadNormal(normals,2,coord,0,1,5,4);

  // Face 1265
  computeFaceQuadNormal(normals,3,coord,1,2,6,5);

  // Face 3762
  computeFaceQuadNormal(normals,4,coord,3,7,6,2);

  // Face 0473
  computeFaceQuadNormal(normals,5,coord,0,4,7,3);

}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeTetraUnitNormals(GeometricEntity& cell,
                                                vector<RealVector>& unitNormals)
{
  /// @todo This information must be gotten from a centralized ElementType data
  const CFuint nbNodesInElem = 4;
  const CFuint nbFacesInElem = 4;

  cf_assert(unitNormals.size() == nbFacesInElem);

  const vector<Node*>* const nodes = cell.getNodes();

  cf_assert(nodes->size() == nbNodesInElem);

  RealMatrix coordinates(nbNodesInElem,DIM_3D);
  RealMatrix faceNormals(nbFacesInElem,DIM_3D);

  putNodesIntoCoord(nodes,coordinates);

  computeTetraNormals(coordinates,faceNormals);

  adimensionalize3d(faceNormals,unitNormals);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::adimensionalize3d(const RealMatrix faceNormals,
                                          vector<RealVector>& unitNormals)
{
  for (CFuint iNorm = 0; iNorm < unitNormals.size(); ++iNorm) {
    const CFreal sum = faceNormals(iNorm,XX) * faceNormals(iNorm,XX) +
                       faceNormals(iNorm,YY) * faceNormals(iNorm,YY) +
                       faceNormals(iNorm,ZZ) * faceNormals(iNorm,ZZ);

    const CFreal invNorm = 1.0 / sqrt(sum);

    unitNormals[iNorm][XX] = faceNormals(iNorm,XX) * invNorm;
    unitNormals[iNorm][YY] = faceNormals(iNorm,YY) * invNorm;
    unitNormals[iNorm][ZZ] = faceNormals(iNorm,ZZ) * invNorm;
  }
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::adimensionalize2d(const RealMatrix faceNormals,
                                          vector<RealVector>& unitNormals)
{
  for (CFuint iNorm = 0; iNorm < unitNormals.size(); ++iNorm) {
    const CFreal sum = faceNormals(iNorm,XX) * faceNormals(iNorm,XX) +
                       faceNormals(iNorm,YY) * faceNormals(iNorm,YY);

    const CFreal invNorm = 1.0 / sqrt(sum);

    unitNormals[iNorm][XX] = faceNormals(iNorm,XX) * invNorm;
    unitNormals[iNorm][YY] = faceNormals(iNorm,YY) * invNorm;
  }
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computePyramUnitNormals(GeometricEntity& cell,
                                                vector<RealVector>& unitNormals)
{
  /// @todo This information must be gotten from a centralized ElementType data
  const CFuint nbNodesInElem = 5;
  const CFuint nbFacesInElem = 5;

  cf_assert(unitNormals.size() == nbFacesInElem);

  const vector<Node*>* const nodes = cell.getNodes();

  cf_assert(nodes->size() == nbNodesInElem);

  RealMatrix coordinates(nbNodesInElem,DIM_3D);
  RealMatrix faceNormals(nbFacesInElem,DIM_3D);

  putNodesIntoCoord(nodes,coordinates);

  computePyramNormals(coordinates,faceNormals);

  adimensionalize3d(faceNormals,unitNormals);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computePrismUnitNormals(GeometricEntity& cell,
                                                vector<RealVector>& unitNormals)
{
  /// @todo This information must be gotten from a centralized ElementType data
  const CFuint nbNodesInElem = 6;
  const CFuint nbFacesInElem = 5;

  cf_assert(unitNormals.size() == nbFacesInElem);

  const vector<Node*>* const nodes = cell.getNodes();

  cf_assert(nodes->size() == nbNodesInElem);

  RealMatrix coordinates(nbNodesInElem,DIM_3D);
  RealMatrix faceNormals(nbFacesInElem,DIM_3D);

  putNodesIntoCoord(nodes,coordinates);

  computePrismNormals(coordinates,faceNormals);

  adimensionalize3d(faceNormals,unitNormals);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeHexaUnitNormals(GeometricEntity& cell,
                                                vector<RealVector>& unitNormals)
{
  /// @todo This information must be gotten from a centralized ElementType data
  const CFuint nbNodesInElem = 8;
  const CFuint nbFacesInElem = 6;

  cf_assert(unitNormals.size() == nbFacesInElem);

  const vector<Node*>* const nodes = cell.getNodes();

  cf_assert(nodes->size() == nbNodesInElem);

  RealMatrix coordinates(nbNodesInElem,DIM_3D);
  RealMatrix faceNormals(nbFacesInElem,DIM_3D);

  putNodesIntoCoord(nodes,coordinates);

  computeHexaNormals(coordinates,faceNormals);

  adimensionalize3d(faceNormals,unitNormals);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeTriagNormals(const RealMatrix& coord,
                                                  RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 3);
  cf_assert(coord.nbCols() == DIM_2D);

  cf_assert(normals.nbRows() == 3);
  cf_assert(normals.nbCols() == DIM_2D);

  // Face 01
  computeFaceLineNormal(normals,0,coord,0,1);

  // Face 12
  computeFaceLineNormal(normals,1,coord,1,2);

  // Face 20
  computeFaceLineNormal(normals,2,coord,2,0);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeTriagUnitNormals(GeometricEntity& cell,
                                                vector<RealVector>& unitNormals)
{
  /// @todo This information must be gotten from a centralized ElementType data
  const CFuint nbNodesInElem = 3;
  const CFuint nbFacesInElem = 3;

  cf_assert(unitNormals.size() == nbFacesInElem);

  const vector<Node*>* const nodes = cell.getNodes();

  cf_assert(nodes->size() == nbNodesInElem);

  RealMatrix coordinates(nbNodesInElem,DIM_2D);
  RealMatrix faceNormals(nbFacesInElem,DIM_2D);

  putNodesIntoCoord(nodes,coordinates);

  computeTriagNormals(coordinates,faceNormals);

  adimensionalize2d(faceNormals,unitNormals);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeQuadNormals(const RealMatrix& coord,
                                                 RealMatrix& normals)
{
  cf_assert(coord.nbRows() == 4);
  cf_assert(coord.nbCols() == DIM_2D);

  cf_assert(normals.nbRows() == 4);
  cf_assert(normals.nbCols() == DIM_2D);

  // Face 01
  computeFaceLineNormal(normals,0,coord,0,1);

  // Face 12
  computeFaceLineNormal(normals,1,coord,1,2);

  // Face 23
  computeFaceLineNormal(normals,2,coord,2,3);

  // Face 30
  computeFaceLineNormal(normals,3,coord,3,0);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::computeQuadUnitNormals(GeometricEntity& cell,
                                               vector<RealVector>& unitNormals)
{
  /// @todo This information must be gotten from a centralized ElementType data
  const CFuint nbNodesInElem = 4;
  const CFuint nbFacesInElem = 4;

  cf_assert(unitNormals.size() == nbFacesInElem);

  const vector<Node*>* const nodes = cell.getNodes();

  cf_assert(nodes->size() == nbNodesInElem);

  RealMatrix coordinates(nbNodesInElem,DIM_2D);
  RealMatrix faceNormals(nbFacesInElem,DIM_2D);

  putNodesIntoCoord(nodes,coordinates);

  computeQuadNormals(coordinates,faceNormals);

  adimensionalize2d(faceNormals,unitNormals);
}

//////////////////////////////////////////////////////////////////////////////

void NormalsCalculator::putNodesIntoCoord(const std::vector<Node*>* const nodes, RealMatrix& coordinates)
{
  for (CFuint n = 0; n < nodes->size(); ++n) {
    coordinates.setRow(*(*nodes)[n],n);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

