// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/BadValueException.hh"

#include "Framework/LocalConnectionData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

LocalConnectionData::LocalConnectionData()
{
  /// Some entries in the tables below are commented out because the corresponding
  /// geometric entities in LocalConnectionDataBuilder are not implemented yet

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::LINE, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofLineOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTriagOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTriagOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER3, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTriagOrder3());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTriagOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTriagOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER3, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTriagOrder3());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER1, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::faceDofTriagOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER2, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::faceDofTriagOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TRIAG, CFPolyOrder::ORDER3, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::faceDofTriagOrder3());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofQuadOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofQuadOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofQuadOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofQuadOrder2());

//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER3, STATE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::faceDofQuadOrder3());
//
//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER4, STATE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::faceDofQuadOrder4());
//
//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER5, STATE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::faceDofQuadOrder5());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER1, STATE, CFPolyForm::SPECTRALFD),
                         LocalConnectionDataBuilder::faceDofQuadOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER2, STATE, CFPolyForm::SPECTRALFD),
                         LocalConnectionDataBuilder::faceDofQuadOrder2());

//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER3, STATE, CFPolyForm::SPECTRALFD),
//                          LocalConnectionDataBuilder::faceDofQuadOrder3());
//
//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER4, STATE, CFPolyForm::SPECTRALFD),
//                          LocalConnectionDataBuilder::faceDofQuadOrder4());
//
//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::QUAD, CFPolyOrder::ORDER5, STATE, CFPolyForm::SPECTRALFD),
//                          LocalConnectionDataBuilder::faceDofQuadOrder5());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTetraOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTetraOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTetraOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofTetraOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER1, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::faceDofTetraOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER2, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::faceDofTetraOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::PYRAM, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofPyramOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::PYRAM, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofPyramOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofPrismOrder1());
  
  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofPrismOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofPrismOrder1());
       
  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofPrismOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofHexaOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::faceDofHexaOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
                         LocalConnectionDataBuilder::faceDofHexaOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
                         LocalConnectionDataBuilder::faceDofHexaOrder2());

//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, NODE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::faceDofHexaOrder3());

//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, STATE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::faceDofHexaOrder3());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
                         LocalConnectionDataBuilder::faceDofHexaOrder2());

//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, NODE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::faceDofHexaOrder3());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER1, STATE, CFPolyForm::SPECTRALFD),
                         LocalConnectionDataBuilder::faceDofHexaOrder1());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, STATE, CFPolyForm::SPECTRALFD),
                         LocalConnectionDataBuilder::faceDofHexaOrder2());

//   _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, STATE, CFPolyForm::SPECTRALFD),
//                          LocalConnectionDataBuilder::faceDofHexaOrder3());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, NODE, CFPolyForm::SERENDIPITY),
                         LocalConnectionDataBuilder::faceDofHexa20NodesOrder2());

  _code2TableFace.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, STATE, CFPolyForm::SERENDIPITY),
                         LocalConnectionDataBuilder::faceDofHexa20NodesOrder2());

  _code2TableFace.sortKeys();


  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofTetraOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofTetraOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofTetraOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofTetraOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER1, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::edgeDofTetraOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::TETRA, CFPolyOrder::ORDER2, STATE, CFPolyForm::SPECTRALFV),
                         LocalConnectionDataBuilder::edgeDofTetraOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::PYRAM, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofPyramOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::PYRAM, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofPyramOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofPrismOrder1());
       
  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofPrismOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofPrismOrder1());
       
  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::PRISM, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofPrismOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofHexaOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, NODE, CFPolyForm::LAGRANGE),
                         LocalConnectionDataBuilder::edgeDofHexaOrder2());

//   _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, NODE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::edgeDofHexaOrder3());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, NODE,CFPolyForm::SERENDIPITY),
                         LocalConnectionDataBuilder::edgeDofHexa20NodesOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER1, STATE, CFPolyForm::LAGRANGE),
       LocalConnectionDataBuilder::edgeDofHexaOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, STATE, CFPolyForm::LAGRANGE),
                         LocalConnectionDataBuilder::edgeDofHexaOrder2());

//   _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, STATE, CFPolyForm::LAGRANGE),
//                          LocalConnectionDataBuilder::edgeDofHexaOrder3());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, STATE, CFPolyForm::SERENDIPITY),
       LocalConnectionDataBuilder::edgeDofHexa20NodesOrder2());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER1, STATE, CFPolyForm::SPECTRALFD),
                             LocalConnectionDataBuilder::edgeDofHexaOrder1());

  _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER2, STATE, CFPolyForm::SPECTRALFD),
                             LocalConnectionDataBuilder::edgeDofHexaOrder2());

//   _code2TableEdge.insert(LocalConnectionData::getCode(CFGeoShape::HEXA, CFPolyOrder::ORDER3, STATE, CFPolyForm::SPECTRALFD),
//                              LocalConnectionDataBuilder::edgeDofHexaOrder3());


  _code2TableEdge.sortKeys();

  prepareMapCode2FaceShape();
  prepareMapCode2EdgeShape();
}

//////////////////////////////////////////////////////////////////////////////

LocalConnectionData::~LocalConnectionData()
{
  // destruct the contents inside the tables
}

//////////////////////////////////////////////////////////////////////////////

LocalConnectionData& LocalConnectionData::getInstance()
{
  static LocalConnectionData theInstance;
  return theInstance;
}

//////////////////////////////////////////////////////////////////////////////

Common::Table<CFuint>* LocalConnectionData::getFaceDofLocal(const CFGeoShape::Type& shape,
                 const CFPolyOrder::Type& geomOrder,
                 const CFDofType& dofType,
                 const CFPolyForm::Type& geomPolyType)
  {
    return _code2TableFace.find(getCode(shape,
            geomOrder,
            dofType,
            geomPolyType));
  }

//////////////////////////////////////////////////////////////////////////////

Common::Table<CFuint>* LocalConnectionData::getEdgeDofLocal(const CFGeoShape::Type& shape,
                 const CFPolyOrder::Type& geomOrder,
                 const CFDofType& dofType,
                 const CFPolyForm::Type& geomPolyType)
  {
    return _code2TableEdge.find(getCode(shape,
            geomOrder,
            dofType,
            geomPolyType));
  }

//////////////////////////////////////////////////////////////////////////////

CFuint LocalConnectionData::getCode(const CFGeoShape::Type& shape,
                        const CFPolyOrder::Type& polyOrder,
                        const CFDofType& dofType,
                        const CFPolyForm::Type& polyType)
  {
    return 10000*(polyOrder+1) + 100*shape + 4*polyType + dofType;
//     return 1000000*(polyOrder+1) + 10000*shape + 100*polyType + dofType;
  }

//////////////////////////////////////////////////////////////////////////////

CFGeoShape::Type LocalConnectionData::getFaceShape(const CFGeoShape::Type& shape,
         const CFuint& iFace)
  {
    return _code2FaceShape.find(getFaceCode(shape, iFace));
  }

//////////////////////////////////////////////////////////////////////////////

CFGeoShape::Type LocalConnectionData::getEdgeShape(const CFGeoShape::Type& shape,
         const CFuint& iEdge)
  {
    return _code2EdgeShape.find(getEdgeCode(shape, iEdge));
  }

//////////////////////////////////////////////////////////////////////////////

CFuint LocalConnectionData::getNbFacesInShape(const CFGeoShape::Type& shape)
{
  switch(shape) {
  case CFGeoShape::LINE:
    return 2;
  case CFGeoShape::TRIAG:
    return 3;
  case CFGeoShape::QUAD:
    return 4;
  case CFGeoShape::TETRA:
    return 4;
  case CFGeoShape::PYRAM:
    return 5;
  case CFGeoShape::PRISM:
    return 5;
  case CFGeoShape::HEXA:
    return 6;
  default:
    std::string msg = std::string("LocalConnectionData::nbFacesInShape() : ") +
      std::string("Shape not defined: ") +
      Common::StringOps::to_str(shape) +
      "\n";
    throw BadValueException(FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint LocalConnectionData::getNbEdgesInShape(const CFGeoShape::Type& shape)
{
  switch(shape) {

  case CFGeoShape::TETRA:
    return 6;
  case CFGeoShape::PYRAM:
    return 8;
  case CFGeoShape::PRISM:
    return 9;
  case CFGeoShape::HEXA:
    return 12;
  default:
    std::string msg = std::string("LocalConnectionData::nbEdgesInShape() : ") +
      std::string("Shape not defined: ") +
      Common::StringOps::to_str(shape) +
      "\n";
    throw BadValueException(FromHere(),msg);
  }
}


//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::setNbFacesPerElement(std::valarray<CFuint>& nbFacesPerElem)
{
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElementTypes = elementType->size();
  CFuint elemID = 0;

  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type curElemType = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID) {
      cf_assert(elemID < nbFacesPerElem.size());
      nbFacesPerElem[elemID] = LocalConnectionData::getNbFacesInShape(curElemType);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::setNbEdgesPerElement(std::valarray<CFuint>& nbEdgesPerElem)
{
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElementTypes = elementType->size();
  CFuint elemID = 0;

  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type curElemType = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID) {
      cf_assert(elemID < nbEdgesPerElem.size());
      nbEdgesPerElem[elemID] = LocalConnectionData::getNbEdgesInShape(curElemType);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::setFaceShapesPerElemType
(vector< vector<CFGeoShape::Type> >& faceShapesPerElemType)
{
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  for (CFuint iElemType = 0; iElemType < nbElemTypes; ++iElemType) {
    const CFuint nbFacesInElem =
      getNbFacesInShape((*elementType)[iElemType].getGeoShape());
    faceShapesPerElemType[iElemType].resize(nbFacesInElem);
  }

  for (CFuint iElemType = 0; iElemType < nbElemTypes; ++iElemType) {
    const CFuint nbFaces = faceShapesPerElemType[iElemType].size();
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      faceShapesPerElemType[iElemType][iFace] =
        getFaceShape((*elementType)[iElemType].getGeoShape(), iFace);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::setEdgeShapesPerElemType
(vector< vector<CFGeoShape::Type> >& edgeShapesPerElemType)
{
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  for (CFuint iElemType = 0; iElemType < nbElemTypes; ++iElemType) {
    const CFuint nbEdgesInElem =
      getNbEdgesInShape((*elementType)[iElemType].getGeoShape());
    edgeShapesPerElemType[iElemType].resize(nbEdgesInElem);
  }

  for (CFuint iElemType = 0; iElemType < nbElemTypes; ++iElemType) {
    const CFuint nbEdges = edgeShapesPerElemType[iElemType].size();
    for (CFuint iEdge = 0; iEdge < nbEdges; ++iEdge) {
      edgeShapesPerElemType[iElemType][iEdge] =
        getEdgeShape((*elementType)[iElemType].getGeoShape(), iEdge);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

CFuint LocalConnectionData::getFaceCode(const CFGeoShape::Type& shape,
       const CFuint& iFace)
  {
    return 10*shape + iFace;
  }

//////////////////////////////////////////////////////////////////////////////

CFuint LocalConnectionData::getEdgeCode(const CFGeoShape::Type& shape,
       const CFuint& iEdge)
  {
    return 10*shape + iEdge;
  }

//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::prepareMapCode2FaceShape()
{
  _code2FaceShape.insert(getFaceCode( CFGeoShape::LINE,0), CFGeoShape::POINT);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::LINE,1), CFGeoShape::POINT);

  _code2FaceShape.insert(getFaceCode( CFGeoShape::TRIAG,0), CFGeoShape::LINE);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::TRIAG,1), CFGeoShape::LINE);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::TRIAG,2), CFGeoShape::LINE);

  _code2FaceShape.insert(getFaceCode( CFGeoShape::QUAD,0), CFGeoShape::LINE);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::QUAD,1), CFGeoShape::LINE);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::QUAD,2), CFGeoShape::LINE);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::QUAD,3), CFGeoShape::LINE);

  _code2FaceShape.insert(getFaceCode( CFGeoShape::TETRA,0), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::TETRA,1), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::TETRA,2), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::TETRA,3), CFGeoShape::TRIAG);

  _code2FaceShape.insert(getFaceCode( CFGeoShape::PYRAM,0), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PYRAM,1), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PYRAM,2), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PYRAM,3), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PYRAM,4), CFGeoShape::TRIAG);

  _code2FaceShape.insert(getFaceCode( CFGeoShape::PRISM,0), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PRISM,1), CFGeoShape::TRIAG);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PRISM,2), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PRISM,3), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::PRISM,4), CFGeoShape::QUAD);

  _code2FaceShape.insert(getFaceCode( CFGeoShape::HEXA,0), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::HEXA,1), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::HEXA,2), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::HEXA,3), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::HEXA,4), CFGeoShape::QUAD);
  _code2FaceShape.insert(getFaceCode( CFGeoShape::HEXA,5), CFGeoShape::QUAD);

  _code2FaceShape.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::prepareMapCode2EdgeShape()
{

  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::TETRA,0), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::TETRA,1), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::TETRA,2), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::TETRA,3), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::TETRA,4), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::TETRA,5), CFGeoShape::LINE);

  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,0), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,1), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,2), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,3), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,4), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,5), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,6), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PYRAM,7), CFGeoShape::LINE);

  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,0), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,1), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,2), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,3), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,4), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,5), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,6), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,7), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::PRISM,8), CFGeoShape::LINE);

  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,0), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,1), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,2), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,3), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,4), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,5), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,6), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,7), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,8), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,9),  CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,10), CFGeoShape::LINE);
  _code2EdgeShape.insert(getEdgeCode( CFGeoShape::HEXA,11), CFGeoShape::LINE);

  _code2EdgeShape.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void LocalConnectionData::print()
{
  CFout << "\n Printing LocalConnectionData::_code2Table :\n";
  _code2TableFace.print();

  CFout << "\n Printing LocalConnectionData::_code2TableEdge :\n";
  _code2TableEdge.print();

  CFout << "\n Printing LocalConnectionData::_code2FaceShape :\n";
  _code2FaceShape.print();

  CFout << "\n Printing LocalConnectionData::_code2EdgeShape :\n";
  _code2EdgeShape.print();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
