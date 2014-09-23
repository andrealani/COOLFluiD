// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Common/Stopwatch.hh"
#include "Common/BadValueException.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "CFmeshCellSplitter/CellSplitter2D.hh"
#include "CFmeshCellSplitter/CFmeshCellSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace CFmeshCellSplitter {

      const CFuint nodesInTriag = 3;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CellSplitter2D,
               MeshFormatConverter,
               CFmeshCellSplitterModule,
               1>
splitter2DProvider("CellSplitter2D");

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

CellSplitter2D::CellSplitter2D(const std::string& name)
  : MeshFormatConverter(name),
    _data(new CFmeshReaderWriterSource()),
    _newTriags(2)
{
   addConfigOptionsTo(this);
  SafePtr<CFmeshReaderWriterSource> ptr = _data.get();
  _reader.setReadData(ptr);
  _writer.setWriteData(ptr);

 }

//////////////////////////////////////////////////////////////////////////////

CellSplitter2D::~CellSplitter2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  for (CFuint i=0;i < _newTriags.size();++i){
    _newTriags[i].resize(3);
  }

}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::checkFormat(const boost::filesystem::path& filepath)
{
//  CFmeshFormatChecker::CFmeshFileChecker checker("checker");
//  checker.check(boost::filesystem::change_extension(filepath,getOriginExtension()));

/// @todo CFmeshFormatChecker should be moved to framework
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::convert(const boost::filesystem::path& fromFilepath,
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
  CFout << "Splitting of the 2D CFmesh took: " << stp.read() << "s\n";
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::convertBack(const boost::filesystem::path& filepath)
{
  CFLog(VERBOSE,"No sence in trying to convert back from a extruded 3D CFmesh" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::split()
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

  transformQuadElements();

  data.setNbElementTypes(1);
  data.consistencyCheck();
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::transformQuadElements()
{
  CFAUTOTRACE;

  convertElementsToTriag();

  migrateElementTypes();
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::convertElementsToTriag()
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

  _pattern.resize(nbTriangles + 2*nbQuads);

  // search for the specific element type IDs
  CFuint elemID = 0;
  for(CFuint iType = 0; iType <  nbElementTypes; ++iType) {

    const CFuint nbElemsPerType = (*elementType)[iType].getNbElems();
    const CFGeoShape::Type elemGeoShape = (*elementType)[iType].getGeoShape();
    for(CFuint iElem = 0; iElem < nbElemsPerType; ++iElem, ++elemID) {

      switch(elemGeoShape) {

      case CFGeoShape::TRIAG:
        _pattern[elemID] = nodesInTriag;
        break;

      case CFGeoShape::QUAD:
        _pattern[elemID] = nodesInTriag;
        _pattern[elemID + nbElemsPerType] = nodesInTriag;
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

  CFuint newNbElements = nbElements;

  elemID = 0;
  CFuint quadElemID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type currShape = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    for(CFuint iElem = 0; iElem < nbElemPerType; ++iElem) {

      switch(currShape) {

      case CFGeoShape::TRIAG:

        for(CFuint localID = 0; localID < 3; ++localID) {

          (*elementNode)(elemID,localID) =
            _oldElemNode (elemID,localID);

          (*elementState)(elemID,localID) =
            _oldElemState(elemID,localID);
        }

        ++elemID;
        break;

      case CFGeoShape::QUAD:
        {
        newNbElements += 1;

        vector<CFuint> newID(2);
        newID[0] = elemID;
        newID[1] = nbElements + quadElemID;

/*        CFuint oldNbElemsPerType = (*elementType)[_triagTypeID].getNbElems();
        (*elementType)[_triagTypeID].setNbElems(oldNbElemsPerType + 1);*/

        vector<CFuint> tempQuad;
        tempQuad.resize(4);
        // Create the Quad

        for(CFuint localID = 0; localID < 4; ++localID) {
          tempQuad[localID] = _oldElemNode(elemID,localID);
        }

        splitQuads(tempQuad);

        //Create the triangles
        for (CFuint iTriag = 0; iTriag < 2; ++iTriag) {
          for(CFuint localID = 0; localID < 3; ++localID) {
            (*elementNode)(newID[iTriag],localID) = _newTriags[iTriag][localID];
            (*elementState)(newID[iTriag],localID) = _newTriags[iTriag][localID];
          }
        }

        ++elemID;
        ++quadElemID;
        }
        break;

      default:

        std::string shape = CFGeoShape::Convert::to_str(currShape);

        std::string msg = std::string("Wrong kind of elements present in 2D mesh: ") +
                       shape +
                       std::string(" ElemID: ") +
                       Common::StringOps::to_str(++elemID);
        throw BadValueException (FromHere(),msg);
        }
    }
  }

  cf_assert(elemID == nbElements);

  data.setNbElements(newNbElements);
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::migrateElementTypes()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr<vector<ElementTypeData> > elementType = data.getElementTypeData();

  // only accept TRIAG and QUAD present
  cf_assert(elementType->size() <= 2);

  elementType->resize(1);
  (*elementType)[0].setGeoShape(CFGeoShape::TRIAG);
  (*elementType)[0].setShape( CFGeoShape::Convert::to_str(CFGeoShape::TRIAG) );
  (*elementType)[0].setNbNodes(nodesInTriag);
  (*elementType)[0].setNbStates(nodesInTriag);
  (*elementType)[0].setNbElems(data.getNbElements());
  (*elementType)[0].setStartIdx(0);
  (*elementType)[0].setGeoOrder(1);
  (*elementType)[0].setSolOrder(1);
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2D::splitQuads(vector<CFuint> quad)
{
  CFAUTOTRACE;

  _newTriags[0][0] = quad[0];
  _newTriags[0][1] = quad[1];
  _newTriags[0][2] = quad[2];

  _newTriags[1][0] = quad[0];
  _newTriags[1][1] = quad[2];
  _newTriags[1][2] = quad[3];

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshCellSplitter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
