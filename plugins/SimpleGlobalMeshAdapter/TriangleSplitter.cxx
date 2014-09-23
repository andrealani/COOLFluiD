// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Common/Stopwatch.hh"
#include "Common/BadValueException.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
// #include "CFmeshFormatChecker/CFmeshFileChecker.hh" // move checker to Framework
#include "SimpleGlobalMeshAdapter/TriangleSplitter.hh"
#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

      const CFuint nodesInTriag = 3;
      const CFuint nodesInQuad = 4;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<TriangleSplitter,
               MeshFormatConverter,
               SimpleGlobalMeshAdapterModule,
               1>
triangleSplitterProvider("TriangleSplitter");

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

TriangleSplitter::TriangleSplitter(const std::string& name)
  : MeshFormatConverter(name),
    _data(new CFmeshReaderWriterSource()),
    _newTriags(3)
{
   addConfigOptionsTo(this);
  SafePtr<CFmeshReaderWriterSource> ptr = _data.get();
  _reader.setReadData(ptr);
  _writer.setWriteData(ptr);

 }

//////////////////////////////////////////////////////////////////////////////

TriangleSplitter::~TriangleSplitter()
{
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  for (CFuint i=0;i < _newTriags.size();++i){
    _newTriags[i].resize(3);
  }

}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::checkFormat(const boost::filesystem::path& filepath)
{
//  CFmeshFormatChecker::CFmeshFileChecker checker("checker");
//  checker.check(boost::filesystem::change_extension(filepath,getOriginExtension()));

/// @todo CFmeshFormatChecker should be moved to framework
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::convert(const boost::filesystem::path& fromFilepath,
       const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  Common::Stopwatch<WallTime> stp;
  stp.start();

  // reads the origin 2D CFmesh
  _reader.readFromFile(fromFilepath);

  // transforms the data
  split();

  // write the new 2D data to the file
  _writer.writeToFile(filepath);

  stp.stop();
  CFout << "Refining of the 2D CFmesh took: " << stp.read() << "s\n";
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::convertBack(const boost::filesystem::path& filepath)
{

  //Nothing to be done here
  cf_assert(false);
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::split()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  cf_assert(data.getGeometricPolyOrder() == CFPolyOrder::ORDER1);
  cf_assert(data.getSolutionPolyOrder()  == CFPolyOrder::ORDER1);

  data.consistencyCheck();

  // set dimension to 2D
  cf_assert(data.getDimension() == DIM_2D);
  data.setDimension(DIM_2D);

  // save connectivities previous to extrusion
  data.copyElementNodeTo(_oldElemNode);
  data.copyElementStateTo(_oldElemState);

  refineElements();

  SafePtr< vector<RealVector> > nodes  = data.getNodeList();
  SafePtr< vector<RealVector> > states = data.getStateList();

  data.setNbUpdatableNodes(nodes->size());
  data.setNbNonUpdatableNodes(0);

  data.setNbUpdatableStates(states->size());
  data.setNbNonUpdatableStates(0);

  data.consistencyCheck();
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::executeRefinement()
{
  CFAUTOTRACE;

  //Mark the elements to refine
  //Mark the refined elements neighbors
  markElements();

  //Split the elements marked for refinement
  refineElements();

  //Split the neighbor elements
  adaptNeighborElements();

  //Split the neighbor elements
  migrateElementTypes();
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::markElements()
{
  CFAUTOTRACE;

///build the cell neighbors

///loop over elements
///Check the cell quality, error or any other...
//  CFmeshReaderWriterSource& data = *(_data.get());

//  const CFuint nbElements = data.getNbElements();

///Mark the cell for refinement and neighbors for adaptation
///(if they are not already marked for refinement themselves)


}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::adaptNeighborElements()
{

  CFAUTOTRACE;

}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::refineElements()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  // this needs to be fixed for FVM
  cf_assert(data.getTotalNbNodes() == data.getTotalNbStates());

  SafePtr< Table<CFuint> > elementNode  = data.getElementNode();
  SafePtr< Table<CFuint> > elementState = data.getElementState();

  elementNode->clear();
  elementState->clear();

  SafePtr< vector<ElementTypeData> > elementType =
    data.getElementTypeData();

  const CFuint nbElements = data.getNbElements();

  CFuint nbQuads = 0;
  CFuint nbTriangles = 0;
  const CFuint nbElementTypes = data.getNbElementTypes();
  cf_assert(nbElementTypes == elementType->size());
  for(CFuint iType = 0; iType <  nbElementTypes; ++iType) {
    if((*elementType)[iType].getGeoShape() ==  CFGeoShape::TRIAG) {
      nbTriangles = (*elementType)[iType].getNbElems();
      _triagTypeID = iType;
    }
    if((*elementType)[iType].getGeoShape() ==  CFGeoShape::QUAD) {
      nbQuads = (*elementType)[iType].getNbElems();
      _quadTypeID = iType;
    }
  }

  ///Compute the number of newTriangles, newQuads
  CFuint nbNewTriangles = 0;
  CFuint nbNewQuads = 0;
  CFuint nbNewNodes = 0;
  CFuint nbNewStates = 0;

  ///Resize the pattern
  _pattern.resize(nbTriangles + nbQuads + nbNewTriangles + nbNewQuads);

  // search for the specific element type IDs
  CFuint elemID = 0;
  for(CFuint iType = 0; iType <  nbElementTypes; ++iType) {

    const CFuint nbElemsPerType = (*elementType)[iType].getNbElems();
    const CFGeoShape::Type elemGeoShape = (*elementType)[iType].getGeoShape();
    for(CFuint iElem = 0; iElem < nbElemsPerType; ++iElem, ++elemID) {

      switch(elemGeoShape) {

      case CFGeoShape::TRIAG:
        _pattern[elemID] = nodesInTriag;
        _pattern[elemID + nbElemsPerType] = nodesInTriag;
        _pattern[elemID + 2*nbElemsPerType] = nodesInTriag;
        break;

      case CFGeoShape::QUAD:
        _pattern[elemID] = nodesInQuad;
        break;

      default:
        std::string shape =
          CFGeoShape::Convert::to_str(elemGeoShape);
        std::string msg = std::string("Wrong kind of elements present in 2D mesh: ") + shape;
        throw BadValueException (FromHere(),msg);
      }
    }
  }

  elementNode->resize(_pattern);
  elementState->resize(_pattern);

  CFuint oldNbNodes = data.getTotalNbNodes();

CFout << "Initial nb Nodes : " <<oldNbNodes <<"\n";
  CFuint nbNodes = oldNbNodes + nbNewNodes;
  CFuint nbStates = oldNbNodes + nbNewStates;

CFout << "Total nb Nodes : " <<nbNodes <<"\n";
std::vector<RealVector> newNodes(nbTriangles);
for(CFuint iNode = 0; iNode < nbTriangles; ++iNode) {
  newNodes[iNode].resize(DIM_2D);
}
  std::vector<CFuint> tempTriag(3);
  elemID = 0;
  CFuint triagElemID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type currShape = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    for(CFuint iElem = 0; iElem < nbElemPerType; ++iElem) {

      switch(currShape) {

      case CFGeoShape::TRIAG:
        {
          newNodes[triagElemID] = 0.;
          for(CFuint localID = 0; localID < 3; ++localID) {
            tempTriag[localID] = _oldElemNode(elemID,localID);
            newNodes[triagElemID] += *(data.getNode(tempTriag[localID]));
          }
          newNodes[triagElemID] /= 3.;
/*CFout << "TriangleNodeID: " <<(tempTriag[0]) <<"\n";
CFout << "TriangleNodeID: " <<(tempTriag[1]) <<"\n";
CFout << "TriangleNodeID: " <<(tempTriag[2]) <<"\n";

CFout << "From: " <<*(data.getNode(tempTriag[0])) <<"\n";
CFout << "From: " <<*(data.getNode(tempTriag[1])) <<"\n";
CFout << "From: " <<*(data.getNode(tempTriag[2])) <<"\n";
CFout << "Created Node: " <<newNodes[triagElemID] <<"\n";*/
          ++elemID;
          ++triagElemID;
        }
        break;

      case CFGeoShape::QUAD:
        {
          ++elemID;
        }
        break;

      default:

        std::string shape = CFGeoShape::Convert::to_str(currShape);

        std::string msg = std::string("Wrong kind of elements present in 2D mesh: ") +
                       shape +
                       std::string(" ElemID: ") +
                       StringOps::to_str(++elemID);
        throw BadValueException (FromHere(),msg);
        }
    }
  }

//Resize nodes and add nodes
//Backup
std::vector<RealVector> oldNodes(oldNbNodes);
for(CFuint iNode = 0; iNode < oldNbNodes; ++iNode) {
  oldNodes[iNode].resize(DIM_2D);
  oldNodes[iNode] = *(data.getNode(iNode));
}

//Resize
data.resizeNodes(nbNodes);
data.resizeStates(nbStates);

//Put back the old nodes
for(CFuint iNode = 0; iNode < oldNbNodes; ++iNode) {
  data.setNode(iNode,oldNodes[iNode]);
}
//Put the new nodes
for(CFuint iNode = 0; iNode < nbTriangles; ++iNode) {
  data.setNode(oldNbNodes + iNode,newNodes[iNode]);
}

data.setNbUpdatableNodes(data.getNbUpdatableNodes() + nbTriangles);
data.setNbUpdatableStates(data.getNbUpdatableStates() + nbTriangles);

CFout << "Initial nbElements: " <<nbElements <<"\n";
  CFuint newNbElements = nbElements;

  elemID = 0;
  triagElemID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type currShape = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    for(CFuint iElem = 0; iElem < nbElemPerType; ++iElem) {

      switch(currShape) {

      case CFGeoShape::TRIAG:
        {

// CFout << "Subdividing triangle:  " <<triagElemID <<"\n";
          newNbElements += 2;

          vector<CFuint> newID(3);
          newID[0] = elemID;
          newID[1] = nbElements + triagElemID;
          newID[2] = nbElements + 2*triagElemID;

  /*        CFuint oldNbElemsPerType = (*elementType)[_triagTypeID].getNbElems();
          (*elementType)[_triagTypeID].setNbElems(oldNbElemsPerType + 1);*/
//   CFout << "Temp triangle\n";
          vector<CFuint> tempTriag;
          tempTriag.resize(3);
          // Create the Triag
          for(CFuint localID = 0; localID < 3; ++localID) {
            tempTriag[localID] = _oldElemNode(elemID,localID);
          }

          CFuint newNodeID = oldNbNodes + triagElemID;
// CFout << "Split triangle\n";
          splitTriags(tempTriag, newNodeID);

          //Create the triangles
          for (CFuint iTriag = 0; iTriag < 3; ++iTriag) {
            for(CFuint localID = 0; localID < 3; ++localID) {
              (*elementNode)(newID[iTriag],localID) = _newTriags[iTriag][localID];
              (*elementState)(newID[iTriag],localID) = _newTriags[iTriag][localID];
            }
          }
//   CFout << "Done with triag\n";
          ++elemID;
          ++triagElemID;
        }
        break;

      case CFGeoShape::QUAD:
        {
          for(CFuint localID = 0; localID < 4; ++localID) {

            (*elementNode)(elemID,localID) =
              _oldElemNode (elemID,localID);

            (*elementState)(elemID,localID) =
              _oldElemState(elemID,localID);
          }
          ++elemID;
        }
        break;

      default:

        std::string shape = CFGeoShape::Convert::to_str(currShape);

        std::string msg = std::string("Wrong kind of elements present in 2D mesh: ") +
                       shape +
                       std::string(" ElemID: ") +
                       StringOps::to_str(++elemID);
        throw BadValueException (FromHere(),msg);
        }
    }
  }

  cf_assert(elemID == nbElements);

  data.setNbElements(newNbElements);
CFout << "New nbElements: " <<newNbElements <<"\n";
}

//////////////////////////////////////////////////////////////////////////////

void TriangleSplitter::migrateElementTypes()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr<vector<ElementTypeData> > elementType = data.getElementTypeData();

  const CFuint nbTypeShapes = elementType->size();

  // only accept TRIAG and QUAD present
  cf_assert(nbTypeShapes <= 2);

  elementType->resize(1);
  (*elementType)[0].setGeoShape(CFGeoShape::TRIAG);
  (*elementType)[0].setShape(CFGeoShape::Convert::to_str(CFGeoShape::TRIAG));
  (*elementType)[0].setNbNodes(nodesInTriag);
  (*elementType)[0].setNbStates(nodesInTriag);
  (*elementType)[0].setNbElems(data.getNbElements());
  (*elementType)[0].setStartIdx(0);
  (*elementType)[0].setGeoOrder(1);
  (*elementType)[0].setSolOrder(1);
}

//////////////////////////////////////////////////////////////////////////////

// void Extruder2D::updateTRSData()
// {
//   CFAUTOTRACE;
//
//   CFuint nbNodesInGeo = 3;
//   CFuint nbStatesInGeo = 3;
//
//   CFmeshReaderWriterSource& data = *(_data.get());
//
//   SafePtr< vector<CFuint> > nbTRs = data.getNbTRs();
//   SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR = data.getNbGeomEntsPerTR();
//
//   for(CFuint iTRS = 0; iTRS < data.getNbTRSs(); ++iTRS ) {
//
//     TRGeoConn& geoConn = data.getTRGeoConn(iTRS);
//
//     for(CFuint iTR = 0; iTR < (*nbTRs)[iTRS]; ++iTR ) {
//
//       GeoConn& geoC = geoConn[iTR];
//       const CFuint oldNbGeos = geoC.size();
//       CFuint nbGeos = oldNbGeos;
//       nbGeos *= 3;
//
//       // make a copy of the old connectivity
//       GeoConn oldConn = geoC;
//
//       geoC.resize(oldNbGeos * 3);
//
//       for(CFuint iLayer = 0; iLayer < _nbLayers; ++iLayer) {
//
//     for (CFuint iGeo = 0; iGeo < oldNbGeos; ++iGeo) {
//
//     // Dirty trick to avoid out of range problems when not splitting
//     CFuint temp = oldNbGeos;
//     if (!_split) temp = 0;
//
//     std::valarray<CFuint>& nodeCon  = geoC[iGeo + (iLayer * nbGeos)].first;
//     std::valarray<CFuint>& stateCon = geoC[iGeo + (iLayer * nbGeos)].second;
//     std::valarray<CFuint>& nodeCon2  = geoC[iGeo+temp + (iLayer * nbGeos)].first;
//     std::valarray<CFuint>& stateCon2 = geoC[iGeo+temp + (iLayer * nbGeos)].second;
//
//     std::valarray<CFuint>& oldNodeCon  = oldConn[iGeo].first;
//     std::valarray<CFuint>& oldStateCon = oldConn[iGeo].second;
//
//     nodeCon.resize(nbNodesInGeo);
//     stateCon.resize(nbStatesInGeo);
//     nodeCon2.resize(nbNodesInGeo);
//     stateCon2.resize(nbStatesInGeo);
//
//     // change node connectivity
//     if (!_split){
//     for(CFuint nodeID = 0; nodeID < halfNodesInGeo; ++nodeID) {
//       nodeCon[nodeID] =
//         oldNodeCon[nodeID] + iLayer * _nbNodesPerLayer;
//       nodeCon[nodeID + halfNodesInGeo] =
//         oldNodeCon[nodeID] + (iLayer + 1) * _nbNodesPerLayer;
//     }
//
//     // change state connectivity
//     for(CFuint stateID = 0; stateID < halfStatesInGeo; ++stateID) {
//       stateCon[stateID] =
//         oldStateCon[stateID] + iLayer * _nbStatesPerLayer;
//       stateCon[stateID + halfStatesInGeo] =
//         oldStateCon[stateID] + (iLayer + 1) * _nbStatesPerLayer;
//     }
//     }
//     else{
//     //CFuint nbNodesInTriag = 3;
//     //CFuint nbNodesInTetra = 4;
//
//     // Compute the tetrahedras
//     vector<CFuint> tempQuad(4);
//     // Create the Quads
//     for(CFuint nodeID = 0; nodeID < halfNodesInGeo; ++nodeID) {
//       tempQuad[nodeID] = oldNodeCon[nodeID] + iLayer * _nbNodesPerLayer;
//       tempQuad[nodeID+halfNodesInGeo] = oldNodeCon[nodeID] + (iLayer + 1) * _nbNodesPerLayer;
//       }
//     swap(tempQuad[2],tempQuad[3]);
//     splitQuads(tempQuad);
//
//     // Set the node connectivity
//     for(CFuint stateID = 0; stateID < 3; ++stateID) {
//       nodeCon[stateID] = _newTriag[0][stateID];
//       nodeCon2[stateID] = _newTriag[1][stateID];
//       }
//
//     // Do the same for state connectivity
//     // Create the Quads
//     for(CFuint nodeID = 0; nodeID < halfNodesInGeo; ++nodeID) {
//       tempQuad[nodeID] = oldStateCon[nodeID] + iLayer * _nbStatesPerLayer;
//       tempQuad[nodeID+halfNodesInGeo] = oldStateCon[nodeID] + (iLayer + 1) * _nbStatesPerLayer;
//       }
//     swap(tempQuad[2],tempQuad[3]);
//     splitQuads(tempQuad);
//
//     // Set the state connectivity
//     for(CFuint stateID = 0; stateID < 3; ++stateID) {
//       stateCon[stateID] = _newTriag[0][stateID];
//       stateCon2[stateID] = _newTriag[1][stateID];
//       }
//
//     }
//
//     if (_split != true){
//     swap(nodeCon[2],nodeCon[3]);
//     swap(stateCon[2],stateCon[3]);
//     }
//   }
//       }
//
//       (*nbGeomEntsPerTR)[iTRS][iTR] = geoC.size();
//     }
//   }
// }
//
// //////////////////////////////////////////////////////////////////////////////
//
 void TriangleSplitter::splitTriags(vector<CFuint> triag, CFuint newNodeID)
 {
   CFAUTOTRACE;

// CFout << "Subdividing triangle:  " <<triag <<"\n";

   _newTriags[0][0] = triag[0];
   _newTriags[0][1] = triag[1];
   _newTriags[0][2] = newNodeID;

// CFout << "Into:  " <<_newTriags[0] <<"\n";
   _newTriags[1][0] = newNodeID;
   _newTriags[1][1] = triag[1];
   _newTriags[1][2] = triag[2];
//  CFout << "Into:  " <<_newTriags[1] <<"\n";
   _newTriags[2][0] = triag[0];
   _newTriags[2][1] = newNodeID;
   _newTriags[2][2] = triag[2];
 // CFout << "Into:  " <<_newTriags[2] <<"\n";
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
