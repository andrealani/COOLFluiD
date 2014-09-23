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

#include "MathTools/MathFunctions.hh"
#include "Environment/ObjectProvider.hh"

#include "MeshTools/MeshTools.hh"
#include "MeshTools/ConcreteQualityCalculator.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConcreteQualityCalculator,QualityCalculator,MeshToolsModule,1>
concreteQualityCalculatorProvider("Concrete");

//////////////////////////////////////////////////////////////////////////////

ConcreteQualityCalculator::ConcreteQualityCalculator(const std::string& name) :
QualityCalculator(name)
{
  _functionMap.insert ( CFGeoShape::POINT,&ConcreteQualityCalculator::calculatePointQuality);
  _functionMap.insert ( CFGeoShape::LINE,&ConcreteQualityCalculator::calculateLineQuality);
  _functionMap.insert ( CFGeoShape::TRIAG,&ConcreteQualityCalculator::calculateTriangleQuality);
  _functionMap.insert ( CFGeoShape::QUAD,&ConcreteQualityCalculator::calculateQuadQuality);
  _functionMap.insert ( CFGeoShape::TETRA,&ConcreteQualityCalculator::calculateTetraQuality);
  _functionMap.insert ( CFGeoShape::PYRAM,&ConcreteQualityCalculator::calculatePyramidQuality);
  _functionMap.insert ( CFGeoShape::PRISM,&ConcreteQualityCalculator::calculatePrismQuality);
  _functionMap.insert ( CFGeoShape::HEXA,&ConcreteQualityCalculator::calculateHexaQuality);

  _functionMap.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

ConcreteQualityCalculator::~ConcreteQualityCalculator()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConcreteQualityCalculator::computeQuality(GeometricEntity* geoEnt)
{
  _geoEnt = geoEnt;

  CALL_MEMBER_FN ( (*this), ( _functionMap.find ( geoEnt->getShape() ) ) ) ();

  return _quality;
}

//////////////////////////////////////////////////////////////////////////////

void ConcreteQualityCalculator::calculateTriangleQuality()
{

CFreal maxEdgeLength;
CFreal edgeLength;
CFreal halfPerimeter;
CFreal area = _geoEnt->computeVolume();

const std::vector<Node*>& nodes = *(_geoEnt->getNodes());
// Edge 0-1
edgeLength = ((*(nodes[0]))[0] - (*(nodes[1]))[0]) * ((*(nodes[0]))[0] - (*(nodes[1]))[0]);
edgeLength += ((*(nodes[0]))[1] - (*(nodes[1]))[1]) * ((*(nodes[0]))[1] - (*(nodes[1]))[1]);

maxEdgeLength = sqrt(edgeLength);
halfPerimeter = maxEdgeLength;

// Edge 0-2
edgeLength = ((*(nodes[0]))[0] - (*(nodes[2]))[0]) * ((*(nodes[0]))[0] - (*(nodes[2]))[0]);
edgeLength += ((*(nodes[0]))[1] - (*(nodes[2]))[1]) * ((*(nodes[0]))[1] - (*(nodes[2]))[1]);

if(maxEdgeLength < sqrt(edgeLength)) maxEdgeLength = sqrt(edgeLength);
halfPerimeter += sqrt(edgeLength);

// Edge 1-2
edgeLength = ((*(nodes[1]))[0] - (*(nodes[2]))[0]) * ((*(nodes[1]))[0] - (*(nodes[2]))[0]);
edgeLength += ((*(nodes[1]))[1] - (*(nodes[2]))[1]) * ((*(nodes[1]))[1] - (*(nodes[2]))[1]);

if(maxEdgeLength < sqrt(edgeLength)) maxEdgeLength = sqrt(edgeLength);
halfPerimeter += sqrt(edgeLength);

halfPerimeter *= 0.5;

_quality = (sqrt(3)/6)*maxEdgeLength*halfPerimeter/area;

}

//////////////////////////////////////////////////////////////////////////////

void ConcreteQualityCalculator::calculateTetraQuality()
{

CFAUTOTRACE;

CFreal maxEdgeLength;
CFreal edgeLength;
CFreal sumFaceAreas;
CFreal volume = _geoEnt->computeVolume();

// Compute the surface of the faces
// faces are not built -> need to be computed
// without going through the geometric entity

const std::vector<Node*>& nodes = *(_geoEnt->getNodes());
_faceNodes.resize(3);

_faceNodes[0] = nodes[0];
_faceNodes[1] = nodes[2];
_faceNodes[2] = nodes[1];

sumFaceAreas = std::abs(_volumeCalc.calculate3DTriagVolume(_faceNodes));

_faceNodes[0] = nodes[0];
_faceNodes[1] = nodes[1];
_faceNodes[2] = nodes[3];

sumFaceAreas += std::abs(_volumeCalc.calculate3DTriagVolume(_faceNodes));

_faceNodes[0] = nodes[1];
_faceNodes[1] = nodes[2];
_faceNodes[2] = nodes[3];

sumFaceAreas += std::abs(_volumeCalc.calculate3DTriagVolume(_faceNodes));

_faceNodes[0] = nodes[0];
_faceNodes[1] = nodes[3];
_faceNodes[2] = nodes[2];

sumFaceAreas += std::abs(_volumeCalc.calculate3DTriagVolume(_faceNodes));

// Edge 0-1
edgeLength = ((*(nodes[0]))[0] - (*(nodes[1]))[0]) * ((*(nodes[0]))[0] - (*(nodes[1]))[0]);
edgeLength += ((*(nodes[0]))[1] - (*(nodes[1]))[1]) * ((*(nodes[0]))[1] - (*(nodes[1]))[1]);
edgeLength += ((*(nodes[0]))[2] - (*(nodes[1]))[2]) * ((*(nodes[0]))[2] - (*(nodes[1]))[2]);
maxEdgeLength = edgeLength;

// Edge 0-2
edgeLength = ((*(nodes[0]))[0] - (*(nodes[2]))[0]) * ((*(nodes[0]))[0] - (*(nodes[2]))[0]);
edgeLength += ((*(nodes[0]))[1] - (*(nodes[2]))[1]) * ((*(nodes[0]))[1] - (*(nodes[2]))[1]);
edgeLength += ((*(nodes[0]))[2] - (*(nodes[2]))[2]) * ((*(nodes[0]))[2] - (*(nodes[2]))[2]);
if(maxEdgeLength < edgeLength) maxEdgeLength = edgeLength;

// Edge 0-3
edgeLength = ((*(nodes[0]))[0] - (*(nodes[3]))[0]) * ((*(nodes[0]))[0] - (*(nodes[3]))[0]);
edgeLength += ((*(nodes[0]))[1] - (*(nodes[3]))[1]) * ((*(nodes[0]))[1] - (*(nodes[3]))[1]);
edgeLength += ((*(nodes[0]))[2] - (*(nodes[3]))[2]) * ((*(nodes[0]))[2] - (*(nodes[3]))[2]);
if(maxEdgeLength < edgeLength) maxEdgeLength = edgeLength;

// Edge 1-2
edgeLength = ((*(nodes[1]))[0] - (*(nodes[2]))[0]) * ((*(nodes[1]))[0] - (*(nodes[2]))[0]);
edgeLength += ((*(nodes[1]))[1] - (*(nodes[2]))[1]) * ((*(nodes[1]))[1] - (*(nodes[2]))[1]);
edgeLength += ((*(nodes[1]))[2] - (*(nodes[2]))[2]) * ((*(nodes[1]))[2] - (*(nodes[2]))[2]);
if(maxEdgeLength < edgeLength) maxEdgeLength = edgeLength;

// Edge 1-3
edgeLength = ((*(nodes[1]))[0] - (*(nodes[3]))[0]) * ((*(nodes[1]))[0] - (*(nodes[3]))[0]);
edgeLength += ((*(nodes[1]))[1] - (*(nodes[3]))[1]) * ((*(nodes[1]))[1] - (*(nodes[3]))[1]);
edgeLength += ((*(nodes[1]))[2] - (*(nodes[3]))[2]) * ((*(nodes[1]))[2] - (*(nodes[3]))[2]);
if(maxEdgeLength < edgeLength) maxEdgeLength = edgeLength;

// Edge 2-3
edgeLength = ((*(nodes[2]))[0] - (*(nodes[3]))[0]) * ((*(nodes[2]))[0] - (*(nodes[3]))[0]);
edgeLength += ((*(nodes[2]))[1] - (*(nodes[3]))[1]) * ((*(nodes[2]))[1] - (*(nodes[3]))[1]);
edgeLength += ((*(nodes[2]))[2] - (*(nodes[3]))[2]) * ((*(nodes[2]))[2] - (*(nodes[3]))[2]);
if(maxEdgeLength < edgeLength) maxEdgeLength = edgeLength;

CFreal radius = 3*volume/sumFaceAreas;
CFreal alpha = sqrt(6)/12;
_quality = alpha*sqrt(maxEdgeLength)/radius;

}

//////////////////////////////////////////////////////////////////////////////

void ConcreteQualityCalculator::calculateQuadQuality()
{

CFAUTOTRACE;

CFreal maxAngle = 0.;
CFreal angleCos;

RealVector edge1(DIM_2D);
RealVector edge2(DIM_2D);

///First check that volume is non-negative
cf_assert(_geoEnt->computeVolume() > 0.);

const std::vector<Node*>& nodes = *(_geoEnt->getNodes());

// Angle between edge 01 and edge 12
edge1 = (*(nodes[1])) - (*(nodes[0]));
edge2 = (*(nodes[2])) - (*(nodes[1]));

edge1.normalize();
edge2.normalize();

angleCos = MathTools::MathFunctions::innerProd(edge1, edge2);
maxAngle = std::max(std::acos(angleCos), maxAngle);

// Angle between edge 12 and edge 23
edge1 = (*(nodes[2])) - (*(nodes[1]));
edge2 = (*(nodes[3])) - (*(nodes[2]));

edge1.normalize();
edge2.normalize();

angleCos = MathTools::MathFunctions::innerProd(edge1, edge2);
maxAngle = std::max(std::acos(angleCos), maxAngle);

// Angle between edge 23 and edge 30
edge1 = (*(nodes[3])) - (*(nodes[2]));
edge2 = (*(nodes[0])) - (*(nodes[3]));

edge1.normalize();
edge2.normalize();

angleCos = MathTools::MathFunctions::innerProd(edge1, edge2);
maxAngle = std::max(std::acos(angleCos), maxAngle);

// Angle between edge 30 and edge 01
edge1 = (*(nodes[0])) - (*(nodes[3]));
edge2 = (*(nodes[1])) - (*(nodes[0]));

edge1.normalize();
edge2.normalize();

angleCos = MathTools::MathFunctions::innerProd(edge1, edge2);
maxAngle = std::max(std::acos(angleCos), maxAngle);


///Define quality from 1 to infinity from the maximum angle.
_quality = 0.5/(1. - (maxAngle/MathTools::MathConsts::CFrealPi()));


}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

