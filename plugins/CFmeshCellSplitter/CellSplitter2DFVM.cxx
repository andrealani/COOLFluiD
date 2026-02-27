// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Common/Stopwatch.hh"
#include "Common/BadValueException.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "CFmeshCellSplitter/CellSplitter2DFVM.hh"
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

Environment::ObjectProvider<CellSplitter2DFVM,
               MeshFormatConverter,
               CFmeshCellSplitterModule,
               1>
splitter2DFVMProvider("CellSplitter2DFVM");

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

CellSplitter2DFVM::CellSplitter2DFVM(const std::string& name)
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

CellSplitter2DFVM::~CellSplitter2DFVM()
{
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  for (CFuint i=0;i < _newTriags.size();++i){
    _newTriags[i].resize(3);
  }

}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::checkFormat(const boost::filesystem::path& filepath)
{
//  CFmeshFormatChecker::CFmeshFileChecker checker("checker");
//  checker.check(boost::filesystem::change_extension(filepath,getOriginExtension()));

/// @todo CFmeshFormatChecker should be moved to framework
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::convert(const boost::filesystem::path& fromFilepath,
                                const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;

  Common::Stopwatch<WallTime> stp;
  stp.start();

  // reads the origin 2D CFmesh
  _reader.readFromFile(fromFilepath);

  // transforms the data
  split();

  if ( fromFilepath.string() == filepath.string() )
  {
#ifdef CF_HAVE_BOOST_1_85
    boost::filesystem::path new_path (  filepath.stem().string() + "_splitted" + filepath.extension().string() );
    boost::filesystem::path new_path2 = filepath.parent_path() / new_path;
#else
  boost::filesystem::path new_path (  basename(filepath) + "_splitted" + extension(filepath) );
  boost::filesystem::path new_path2 = filepath.branch_path() / new_path;
#endif
  }

  // write the new 2D data to the file
  _writer.writeToFile(filepath);

  stp.stop();
  CFout << "Splitting of the 2D CFmesh took: " << stp.read() << "s\n";
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::convertBack(const boost::filesystem::path& filepath)
{
  CFLog(VERBOSE,"No sence in trying to convert back from a extruded 3D CFmesh" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::split()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  cf_assert(data.getGeometricPolyOrder() == CFPolyOrder::ORDER1);
//  cf_assert(data.getSolutionPolyOrder()  == CFPolyOrder::ORDER0);

  data.consistencyCheck();

  // set dimension to 2D
  cf_assert(data.getDimension() == DIM_2D);
  data.setDimension(DIM_2D);

  // save connectivities previous to extrusion
  data.copyElementNodeTo(_oldElemNode);
  data.copyElementStateTo(_oldElemState);

  transformQuadElements();
  
  SafePtr< vector<CFreal> > nodes  = data.getNodeList();
  SafePtr< vector<CFreal> > states = data.getStateList();
  
  data.setNbUpdatableNodes(nodes->size()/data.getDimension());
  data.setNbNonUpdatableNodes(0);
  
  data.setNbUpdatableStates(_statePattern.size());
  data.setNbNonUpdatableStates(0);
  
  data.setNbElementTypes(1);
  data.consistencyCheck();
  data.setWithSolution(false);
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::transformQuadElements()
{
  CFAUTOTRACE;

  convertElementsToTriag();

  convertTRSConnectivity();

  migrateElementTypes();
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::convertElementsToTriag()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

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

  _nodePattern.resize(nbTriangles + 2*nbQuads);
  _statePattern.resize(nbTriangles + 2*nbQuads);

  // search for the specific element type IDs
  CFuint elemID = 0;
  CFuint quadElemID = 0;
  for(CFuint iType = 0; iType <  nbElementTypes; ++iType) {

    const CFuint nbElemsPerType = (*elementType)[iType].getNbElems();
    const CFGeoShape::Type elemGeoShape = (*elementType)[iType].getGeoShape();
    for(CFuint iElem = 0; iElem < nbElemsPerType; ++iElem, ++elemID) {

      switch(elemGeoShape) {

      case CFGeoShape::TRIAG:
        _nodePattern[elemID] = nodesInTriag;
        _statePattern[elemID] = 1;
        break;

      case CFGeoShape::QUAD:
        _nodePattern[elemID] = nodesInTriag;
        _nodePattern[elemID + nbElemsPerType] = nodesInTriag;
        _statePattern[elemID] = 1;
        _statePattern[nbElements + quadElemID] = 1;
        ++quadElemID;
        break;

      default:
        std::string shape =
          CFGeoShape::Convert::to_str(elemGeoShape);
        std::string msg = std::string("Wrong kind of elements present in 2D mesh: ") + shape;
        throw BadValueException (FromHere(),msg);
      }
    }
  }

  elementNode->resize(_nodePattern);
  elementState->resize(_statePattern);

  CFuint newNbElements = nbElements;

  elemID = 0;
  quadElemID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type currShape = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    for(CFuint iElem = 0; iElem < nbElemPerType; ++iElem) {

      switch(currShape) {

      case CFGeoShape::TRIAG:

        for(CFuint localID = 0; localID < 3; ++localID) {

          (*elementNode)(elemID,localID) =
            _oldElemNode (elemID,localID);
        }

        // state
        (*elementState)(elemID,0) = elemID;

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
          }
        }

        // state
        (*elementState)(newID[0],0) = elemID;
        (*elementState)(newID[1],0) = nbElements + quadElemID;

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

void CellSplitter2DFVM::migrateElementTypes()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr<vector<ElementTypeData> > elementType = data.getElementTypeData();
  
  // only accept TRIAG and QUAD present
  cf_assert(elementType->size() <= 2);

  elementType->resize(1);
  (*elementType)[0].setGeoShape(CFGeoShape::TRIAG);
  (*elementType)[0].setShape(CFGeoShape::Convert::to_str(CFGeoShape::TRIAG));
  (*elementType)[0].setNbNodes(nodesInTriag);
  (*elementType)[0].setNbStates(1);
  (*elementType)[0].setNbElems(data.getNbElements());
  (*elementType)[0].setStartIdx(0);
  (*elementType)[0].setGeoOrder(1);
  (*elementType)[0].setSolOrder(0);
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::convertTRSConnectivity()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr< vector<ElementTypeData> > elementType =
    data.getElementTypeData();
  SafePtr< Table<CFuint> > elementNode  = data.getElementNode();
  SafePtr< vector<CFuint> > nbTRs = data.getNbTRs();
  SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR = data.getNbGeomEntsPerTR();
  const CFuint nbQuads = (*elementType)[_quadTypeID].getNbElems();

  for(CFuint iTRS = 0; iTRS < data.getNbTRSs(); ++iTRS ) {

    TRGeoConn& geoConn = data.getTRGeoConn(iTRS);

    for(CFuint iTR = 0; iTR < (*nbTRs)[iTRS]; ++iTR ) {

      GeoConn& geoC = geoConn[iTR];
      const CFuint nbGeos = geoC.size();

      // make a copy of the old connectivity
      GeoConn oldConn = geoC;

      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {

        std::valarray<CFuint>& nodeCon = geoC[iGeo].first;
        std::valarray<CFuint>& stateCon = geoC[iGeo].second;

        //look element with state
        CFuint match=0;
        for (CFuint iElemNode = 0; iElemNode < 3; ++iElemNode) {
          for (CFuint iFaceNode = 0; iFaceNode < 2; ++iFaceNode) {
            if((*elementNode)(stateCon[0],iElemNode) == nodeCon[iFaceNode]) match++;
          }
        }

        //Modify state connectivity
        if(match != 2) stateCon[0] += nbQuads;

        //check that now it is correct
        match=0;
        for (CFuint iElemNode = 0; iElemNode < 3; ++iElemNode) {
          for (CFuint iFaceNode = 0; iFaceNode < 2; ++iFaceNode) {
            if((*elementNode)(stateCon[0],iElemNode) == nodeCon[iFaceNode]) match++;
          }
        }
        cf_assert(match==2);

      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter2DFVM::splitQuads(vector<CFuint> quad)
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
