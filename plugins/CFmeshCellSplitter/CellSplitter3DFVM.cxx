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
#include "CFmeshCellSplitter/CellSplitter3DFVM.hh"
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
      const CFuint nodesInQuad  = 4;

      const CFuint nodesInTetra = 4;
      const CFuint nodesInPyram = 5;
      const CFuint nodesInPrism = 6;
      const CFuint nodesInHexa  = 8;

      const CFuint statesInTetra = 1;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CellSplitter3DFVM,
               MeshFormatConverter,
               CFmeshCellSplitterModule,
               1>
splitter3DFVMProvider("CellSplitter3DFVM");

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

CellSplitter3DFVM::CellSplitter3DFVM(const std::string& name)
  : MeshFormatConverter(name),
    _data(new CFmeshReaderWriterSource()),
    _newTetras(6),
    _newTriags(2)
{
   addConfigOptionsTo(this);
  SafePtr<CFmeshReaderWriterSource> ptr = _data.get();
  _reader.setReadData(ptr);
  _writer.setWriteData(ptr);

 }

//////////////////////////////////////////////////////////////////////////////

CellSplitter3DFVM::~CellSplitter3DFVM()
{
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  for (CFuint i=0; i < _newTetras.size(); ++i){
    _newTetras[i].resize(4);
  }

  for (CFuint i=0;i < _newTriags.size();++i){
    _newTriags[i].resize(3);
  }

}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::checkFormat(const boost::filesystem::path& filepath)
{
//  CFmeshFormatChecker::CFmeshFileChecker checker("checker");
//  checker.check(boost::filesystem::change_extension(filepath,getOriginExtension()));

/// @todo CFmeshFormatChecker should be moved to framework
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::convert(const boost::filesystem::path& fromFilepath,
       const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  Common::Stopwatch<WallTime> stp;
  stp.start();

  // reads the original 3D CFmesh
  _reader.readFromFile(fromFilepath);

  // transforms the data
  split();

//   if ( fromFilepath.string() == filepath.string() ) {
//     boost::filesystem::path new_path (  basename(filepath) + "_splitted" + extension(filepath) );
//     boost::filesystem::path new_path2 = filepath.branch_path() / new_path;
//   }

  // write the new 3D data to the file
  _writer.writeToFile(filepath);

  stp.stop();
  CFout << "Splitting of the 3D CFmesh took: " << stp.read() << "s\n";
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::convertBack(const boost::filesystem::path& filepath)
{
  CFLog(VERBOSE,"No chance to convert back (splitting instructions already lost)" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::split()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  // check data properties
  cf_assert(data.getGeometricPolyOrder() == CFPolyOrder::ORDER1);
  cf_assert(data.getSolutionPolyOrder()  == CFPolyOrder::ORDER0);
  cf_assert(data.getDimension() == DIM_3D);

  data.consistencyCheck();

  // save connectivities previous to extrusion
  data.copyElementNodeTo(_oldElemNode);
  data.copyElementStateTo(_oldElemState);

  convertElementsToTetra();

  migrateElementTypes();

  updateTRSData();

  data.setNbElementTypes(1);
  data.setNbUpdatableStates(data.getNbElements());
  data.setNbNonUpdatableStates(0);

  data.consistencyCheck();
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::calculateNbNewElements()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr< vector<ElementTypeData> > elementType  =  data.getElementTypeData();

  // number of all elements as they are (not splitted)
  const CFuint oldNbElements = data.getNbElements();

  // number of types contained in CFmesh (Tetra+Hexa => nbElementTypes=2)
  const CFuint nbElementTypes = data.getNbElementTypes();

  // reads connectivity
  SafePtr< Table<CFuint> > elementNode  = data.getElementNode();
  SafePtr< Table<CFuint> > elementState = data.getElementState();

  CFuint newNbElements = oldNbElements;

  CFuint prevNbElem = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type currShape = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    switch( currShape ) {
      case CFGeoShape::TETRA : // no new elements
                     break;

      case CFGeoShape::PYRAM : // always splitted into 2 tetras (i.e. one new element)
                     newNbElements += nbElemPerType;
                     break;

      case CFGeoShape::PRISM : // always splitted into 3 tetras (i.e. two new elements)
                     newNbElements += 2*nbElemPerType;
                     break;

      case CFGeoShape::HEXA  : // splitted into either 5 or 6 tetras
                     for(CFuint iHexa = 0; iHexa < nbElemPerType; ++iHexa){
                       CFuint newTetraFromHexa = nbNewTetras(iHexa + prevNbElem);
                       newNbElements += ( newTetraFromHexa - 1 );
                     }
                     break;

      default      : std::string shape = CFGeoShape::Convert::to_str(currShape);
                     std::string msg = std::string("Unhandled type of element present in 3D mesh: ") + shape;
                     throw BadValueException (FromHere(),msg);
    }

    prevNbElem += nbElemPerType;
  }
  // resize connectivity list
  data.setNbElements(newNbElements);
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::convertElementsToTetra()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr< vector<ElementTypeData> > elementType  =  data.getElementTypeData();

  const CFuint oldNbElements = data.getNbElements();
  calculateNbNewElements();
  const CFuint newNbElements = data.getNbElements();

  SafePtr< Table<CFuint> > elementNode  = data.getElementNode();
  SafePtr< Table<CFuint> > elementState = data.getElementState();

  elementNode->clear();
  elementState->clear();

  const CFuint nbElementTypes = data.getNbElementTypes();

  // search for the specific element type IDs
  for(CFuint iType = 0; iType <  nbElementTypes; ++iType) {

    if((*elementType)[iType].getGeoShape() ==  CFGeoShape::TETRA) { // 4 pts per element
      _tetraTypeID = iType;
    }

    if((*elementType)[iType].getGeoShape() ==  CFGeoShape::PYRAM) { // 5 pts per element
      _pyramTypeID = iType;
    }

    if((*elementType)[iType].getGeoShape() ==  CFGeoShape::PRISM) { // 6 pts per element
      _prismTypeID = iType;
    }

    if((*elementType)[iType].getGeoShape() ==  CFGeoShape::HEXA) {  // 8 pts per element
      _hexaTypeID = iType;
    }
  }

  elementNode->resize(newNbElements, nodesInTetra );
  elementState->resize(newNbElements, statesInTetra);

  CFuint elemID = 0;
  CFuint shiftID = 0;

  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFGeoShape::Type currShape = (*elementType)[iType].getGeoShape();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    vector<CFuint> newNodeID;
    vector<CFuint> tempElement;
    for(CFuint iElem = 0; iElem < nbElemPerType; ++iElem) {

      switch(currShape) {

        case CFGeoShape::TETRA:
          cf_assert(iType == _tetraTypeID);
          (*elementState)(elemID,0) = elemID;
          for(CFuint localID = 0; localID < 4; ++localID) {
            // no change - copy old connectivity
            (*elementNode)(elemID,localID) = _oldElemNode( elemID, localID );
          }
          ++elemID;
          break;

        case CFGeoShape::PYRAM:
          cf_assert(iType == _pyramTypeID);

          newNodeID.resize(_newTetrasSize);

          newNodeID[0] = elemID;
          newNodeID[1] = oldNbElements + shiftID;

          tempElement.resize(5);

          for(CFuint localID=0; localID < 5; ++localID){
            tempElement[localID] = _oldElemNode( elemID, localID );
          }

          splitPyram(tempElement);

          //Create the tetrahedras
          for (CFuint iTetra = 0; iTetra < _newTetrasSize; ++iTetra) {
            (*elementState)(newNodeID[iTetra],0) = newNodeID[iTetra];

            for(CFuint localID = 0; localID < 4; ++localID) {
              (*elementNode)( newNodeID[iTetra],localID) = _newTetras[iTetra][localID];
            }
          }
          ++shiftID;
          ++elemID;
          break;

        case CFGeoShape::PRISM:
          cf_assert(iType == _prismTypeID);

          newNodeID.resize(_newTetrasSize);

          newNodeID[0] = elemID;
          newNodeID[1] = oldNbElements + shiftID;
          newNodeID[2] = oldNbElements + shiftID + 1;

          tempElement.resize(6);

          for(CFuint localID=0; localID < 6; ++localID){
            tempElement[localID] = _oldElemNode( elemID, localID );
          }

          splitPrism(tempElement);

          //Create the tetrahedras
          for (CFuint iTetra = 0; iTetra < _newTetrasSize; ++iTetra) {
            (*elementState)(newNodeID[iTetra],0) = newNodeID[iTetra];
            for(CFuint localID = 0; localID < 4; ++localID) {
              (*elementNode)( newNodeID[iTetra],localID) = _newTetras[iTetra][localID];
            }
          }
          shiftID += 2 ;
          ++elemID;
          break;

        case CFGeoShape::HEXA:
          cf_assert(iType == _hexaTypeID);

          tempElement.resize(8);

          for(CFuint localID=0; localID < 8; ++localID){
            tempElement[localID] = _oldElemNode( elemID, localID );
          }

          splitHexa(tempElement);

          newNodeID.resize(_newTetrasSize); // how many elements will be generated from the current one

          newNodeID[0] = elemID;
          for(CFuint iTetra = 1; iTetra < _newTetrasSize; ++iTetra) {
            newNodeID[iTetra] = oldNbElements + shiftID + iTetra - 1;
          }

          //Create the tetrahedras
          for (CFuint iTetra = 0; iTetra < _newTetrasSize; ++iTetra) {
            (*elementState)(newNodeID[iTetra],0) = newNodeID[iTetra];
            for(CFuint localID = 0; localID < 4; ++localID) {
              (*elementNode)( newNodeID[iTetra],localID) = _newTetras[iTetra][localID];
            }
          }
          shiftID += (_newTetrasSize - 1) ;
          ++elemID;
          break;

        default:
          std::string shape = CFGeoShape::Convert::to_str(currShape);
          std::string msg = std::string("Wrong kind of elements present in 3D mesh: ") +
                         shape +
                         std::string(" ElemID: ") +
                         Common::StringOps::to_str(++elemID);
          throw BadValueException (FromHere(),msg);
        }
    }
  }

  cf_assert((elemID+shiftID) == newNbElements);

}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::migrateElementTypes()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr<vector<ElementTypeData> > elementType = data.getElementTypeData();

  elementType->resize(1);
  (*elementType)[0].setGeoShape(CFGeoShape::TETRA);
  (*elementType)[0].setShape("Tetra");
  (*elementType)[0].setNbNodes(nodesInTetra);
  (*elementType)[0].setNbStates(statesInTetra);
  (*elementType)[0].setNbElems(data.getNbElements());
  (*elementType)[0].setStartIdx(0);
  (*elementType)[0].setGeoOrder(1);
  (*elementType)[0].setSolOrder(0);
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::splitPyram(vector<CFuint> pyram)
{
  _newTetrasSize = 2;

  // Pyramid in the article notation (that's why '+1')
  vector<CFuint> NodeGlobalID( nodesInPyram + 1 );

  // assign node numbering of the pyramid
  for (CFuint i=0; i<nodesInPyram ; ++i){
    NodeGlobalID[i+1] = pyram[i];
  }

  if ( min( NodeGlobalID[1], NodeGlobalID[3] ) < min( NodeGlobalID[2], NodeGlobalID[4] ) ) {

    _newTetras[0][0] = NodeGlobalID[ 1 ];
    _newTetras[0][1] = NodeGlobalID[ 2 ];
    _newTetras[0][2] = NodeGlobalID[ 3 ];
    _newTetras[0][3] = NodeGlobalID[ 5 ];

    _newTetras[1][0] = NodeGlobalID[ 1 ];
    _newTetras[1][1] = NodeGlobalID[ 3 ];
    _newTetras[1][2] = NodeGlobalID[ 4 ];
    _newTetras[1][3] = NodeGlobalID[ 5 ];
  }

  else{
    if ( min( NodeGlobalID[2], NodeGlobalID[4] ) < min( NodeGlobalID[1], NodeGlobalID[3] ) ) {

      _newTetras[0][0] = NodeGlobalID[ 2 ];
      _newTetras[0][1] = NodeGlobalID[ 3 ];
      _newTetras[0][2] = NodeGlobalID[ 4 ];
      _newTetras[0][3] = NodeGlobalID[ 5 ];

      _newTetras[1][0] = NodeGlobalID[ 2 ];
      _newTetras[1][1] = NodeGlobalID[ 4 ];
      _newTetras[1][2] = NodeGlobalID[ 1 ];
      _newTetras[1][3] = NodeGlobalID[ 5 ];
    }
    else {
      CFout << "Problem in PYRAM Splitting" << "\n";
      for (CFuint i=0; i < 5 ; ++i) CFout << "Pyram: " << pyram[i] << "\n";
      cf_assert(0);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::splitPrism(vector<CFuint> prism)
{
  _newTetrasSize = 3;

  // Prism in the article notation (that's why '+1')
  vector<CFuint> NodeGlobalID( nodesInPrism + 1 );
  // Indirections
  vector<CFuint> ii( nodesInPrism + 1 );

  // assign node numbering of the prism
  for (CFuint i=0; i<nodesInPrism ; ++i){
    NodeGlobalID[i+1] = prism[i];
  }

  _SmallestFaceID = -1;

  getSmallestID(nodesInPrism, NodeGlobalID);

  cf_assert( _SmallestFaceID >= 0 );

  //create indirections
  switch( _SmallestFaceID ) {
    case 1 : ii[1] = 1;
             ii[2] = 2;
             ii[3] = 3;
             ii[4] = 4;
             ii[5] = 5;
             ii[6] = 6;
             break;

    case 2 : ii[1] = 2;
             ii[2] = 3;
             ii[3] = 1;
             ii[4] = 5;
             ii[5] = 6;
             ii[6] = 4;
             break;

    case 3 : ii[1] = 3;
             ii[2] = 1;
             ii[3] = 2;
             ii[4] = 6;
             ii[5] = 4;
             ii[6] = 5;
             break;

    case 4 : ii[1] = 4;
             ii[2] = 6;
             ii[3] = 5;
             ii[4] = 1;
             ii[5] = 3;
             ii[6] = 2;
             break;

    case 5 : ii[1] = 5;
             ii[2] = 4;
             ii[3] = 6;
             ii[4] = 2;
             ii[5] = 1;
             ii[6] = 3;
             break;

    case 6 : ii[1] = 6;
             ii[2] = 5;
             ii[3] = 4;
             ii[4] = 3;
             ii[5] = 2;
             ii[6] = 1;
             break;

    default: CFerr << " When splitting PRISM into TETRA - face index out of range.\n";
             break;
  }

  if ( min( NodeGlobalID[ii[2]], NodeGlobalID[ii[6]] ) < min( NodeGlobalID[ii[3]], NodeGlobalID[ii[5]] ) ) {

    _newTetras[0][0] = NodeGlobalID[ ii[1] ];
    _newTetras[0][1] = NodeGlobalID[ ii[2] ];
    _newTetras[0][2] = NodeGlobalID[ ii[3] ];
    _newTetras[0][3] = NodeGlobalID[ ii[6] ];

    _newTetras[1][0] = NodeGlobalID[ ii[1] ];
    _newTetras[1][1] = NodeGlobalID[ ii[2] ];
    _newTetras[1][2] = NodeGlobalID[ ii[6] ];
    _newTetras[1][3] = NodeGlobalID[ ii[5] ];

    _newTetras[2][0] = NodeGlobalID[ ii[1] ];
    _newTetras[2][1] = NodeGlobalID[ ii[5] ];
    _newTetras[2][2] = NodeGlobalID[ ii[6] ];
    _newTetras[2][3] = NodeGlobalID[ ii[4] ];
  }

  else{
    if ( min( NodeGlobalID[ii[3]],NodeGlobalID[ii[5]] ) < min( NodeGlobalID[ii[2]],NodeGlobalID[ii[6]] ) ) {

      _newTetras[0][0] = NodeGlobalID[ ii[1] ];
      _newTetras[0][1] = NodeGlobalID[ ii[2] ];
      _newTetras[0][2] = NodeGlobalID[ ii[3] ];
      _newTetras[0][3] = NodeGlobalID[ ii[5] ];

      _newTetras[1][0] = NodeGlobalID[ ii[1] ];
      _newTetras[1][1] = NodeGlobalID[ ii[5] ];
      _newTetras[1][2] = NodeGlobalID[ ii[3] ];
      _newTetras[1][3] = NodeGlobalID[ ii[6] ];

      _newTetras[2][0] = NodeGlobalID[ ii[1] ];
      _newTetras[2][1] = NodeGlobalID[ ii[5] ];
      _newTetras[2][2] = NodeGlobalID[ ii[6] ];
      _newTetras[2][3] = NodeGlobalID[ ii[4] ];
    }
    else {
      CFout << "Problem in PRISM Splitting" << "\n";
      for (CFuint i=0; i < 6 ; ++i) CFout << "Prism: " << prism[i] << "\n";
      cf_assert(0);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::splitHexa(vector<CFuint> hexa)
{
  // Hexa in the article notation (that's why '+1')
  vector<CFuint> NodeGlobalID( nodesInHexa + 1 );
  // Indirections
  vector<CFuint> ii( nodesInHexa + 1 );

  // assign node numbering of the hexahedron
  for (CFuint i=0; i<nodesInHexa ; ++i){
    NodeGlobalID[i+1] = hexa[i];
  }

  _SmallestFaceID = -1;
  getSmallestID(nodesInHexa, NodeGlobalID);
  cf_assert( _SmallestFaceID >= 0 );

  // create indirections - rotate hexa in order to have smallest index in lower-front-left corner
  switch( _SmallestFaceID ) {
    case 1 : ii[1] = 1;
             ii[2] = 2;
             ii[3] = 3;
             ii[4] = 4;
             ii[5] = 5;
             ii[6] = 6;
             ii[7] = 7;
             ii[8] = 8;
             break;

    case 2 : ii[1] = 2;
             ii[2] = 1;
             ii[3] = 5;
             ii[4] = 6;
             ii[5] = 3;
             ii[6] = 4;
             ii[7] = 8;
             ii[8] = 7;
             break;

    case 3 : ii[1] = 3;
             ii[2] = 2;
             ii[3] = 6;
             ii[4] = 7;
             ii[5] = 4;
             ii[6] = 1;
             ii[7] = 5;
             ii[8] = 8;
             break;

    case 4 : ii[1] = 4;
             ii[2] = 1;
             ii[3] = 2;
             ii[4] = 3;
             ii[5] = 8;
             ii[6] = 5;
             ii[7] = 6;
             ii[8] = 7;
             break;

    case 5 : ii[1] = 5;
             ii[2] = 1;
             ii[3] = 4;
             ii[4] = 8;
             ii[5] = 6;
             ii[6] = 2;
             ii[7] = 3;
             ii[8] = 7;
             break;

    case 6 : ii[1] = 6;
             ii[2] = 2;
             ii[3] = 1;
             ii[4] = 5;
             ii[5] = 7;
             ii[6] = 3;
             ii[7] = 4;
             ii[8] = 8;
             break;

    case 7 : ii[1] = 7;
             ii[2] = 3;
             ii[3] = 2;
             ii[4] = 6;
             ii[5] = 8;
             ii[6] = 4;
             ii[7] = 1;
             ii[8] = 5;
             break;

    case 8 : ii[1] = 8;
             ii[2] = 4;
             ii[3] = 3;
             ii[4] = 7;
             ii[5] = 5;
             ii[6] = 1;
             ii[7] = 2;
             ii[8] = 6;
             break;

    default: CFerr << " When splitting HEXA into TETRA - face index out of range.\n";
             break;
  }

  // compute how many diagonals pass through opposite corner
  CFuint DiagCount = 0;
  if ( min( NodeGlobalID[ ii[2] ], NodeGlobalID[ ii[7] ] ) < min( NodeGlobalID[ ii[3] ], NodeGlobalID[ ii[6] ]) )
    DiagCount += 100;
  if ( min( NodeGlobalID[ ii[4] ], NodeGlobalID[ ii[7] ] ) < min( NodeGlobalID[ ii[3] ], NodeGlobalID[ ii[8] ]) )
    DiagCount += 10;
  if ( min( NodeGlobalID[ ii[5] ], NodeGlobalID[ ii[7] ] ) < min( NodeGlobalID[ ii[6] ], NodeGlobalID[ ii[8] ]) )
    DiagCount += 1;

  CFuint nbDiagonals;
  CFuint tmpII;
  // explicit number of diagonals
  if( DiagCount  == 0 )                                                  nbDiagonals = 0;
  if( (DiagCount == 1 ) || (DiagCount == 10)  || (DiagCount == 100) )    nbDiagonals = 1;
  if( (DiagCount == 11) || (DiagCount == 101) || (DiagCount == 110) )    nbDiagonals = 2;
  if( DiagCount  == 111 )                                                nbDiagonals = 3;

  // rotating hexahedron along  <-Corner1--Corner7->  to have suitable face numbering
  if( (DiagCount == 0) || (DiagCount == 11) || (DiagCount == 100) || (DiagCount == 111) ){
    // no rotation needed
  }
  if( (DiagCount == 1) || (DiagCount == 110) ){
    // rotate 120 CounterClockWise
    tmpII = ii[2];
    ii[2] = ii[5];
    ii[5] = ii[4];
    ii[4] = tmpII;
    tmpII = ii[6];
    ii[6] = ii[8];
    ii[8] = ii[3];
    ii[3] = tmpII;
  }
  if( (DiagCount == 10) || (DiagCount == 101) ){
    // rotate 240 CounterClockWise
    tmpII = ii[2];
    ii[2] = ii[4];
    ii[4] = ii[5];
    ii[5] = tmpII;
    tmpII = ii[6];
    ii[6] = ii[3];
    ii[3] = ii[8];
    ii[8] = tmpII;
  }

  switch( nbDiagonals ){
    case 0 : // no diagonals passing through vertex 7
             _newTetras[0][0] = NodeGlobalID[ ii[1] ];
             _newTetras[0][1] = NodeGlobalID[ ii[2] ];
             _newTetras[0][2] = NodeGlobalID[ ii[3] ];
             _newTetras[0][3] = NodeGlobalID[ ii[6] ];

             _newTetras[1][0] = NodeGlobalID[ ii[1] ];
             _newTetras[1][1] = NodeGlobalID[ ii[3] ];
             _newTetras[1][2] = NodeGlobalID[ ii[8] ];
             _newTetras[1][3] = NodeGlobalID[ ii[6] ];

             _newTetras[2][0] = NodeGlobalID[ ii[1] ];
             _newTetras[2][1] = NodeGlobalID[ ii[3] ];
             _newTetras[2][2] = NodeGlobalID[ ii[4] ];
             _newTetras[2][3] = NodeGlobalID[ ii[8] ];

             _newTetras[3][0] = NodeGlobalID[ ii[1] ];
             _newTetras[3][1] = NodeGlobalID[ ii[6] ];
             _newTetras[3][2] = NodeGlobalID[ ii[8] ];
             _newTetras[3][3] = NodeGlobalID[ ii[5] ];

             _newTetras[4][0] = NodeGlobalID[ ii[3] ];
             _newTetras[4][1] = NodeGlobalID[ ii[8] ];
             _newTetras[4][2] = NodeGlobalID[ ii[6] ];
             _newTetras[4][3] = NodeGlobalID[ ii[7] ];

             _newTetrasSize = 5;
            break;

    case 1 : // 1 diagonal passing through vertex 7
             _newTetras[0][0] = NodeGlobalID[ ii[1] ];
             _newTetras[0][1] = NodeGlobalID[ ii[6] ];
             _newTetras[0][2] = NodeGlobalID[ ii[8] ];
             _newTetras[0][3] = NodeGlobalID[ ii[5] ];

             _newTetras[1][0] = NodeGlobalID[ ii[1] ];
             _newTetras[1][1] = NodeGlobalID[ ii[2] ];
             _newTetras[1][2] = NodeGlobalID[ ii[8] ];
             _newTetras[1][3] = NodeGlobalID[ ii[6] ];

             _newTetras[2][0] = NodeGlobalID[ ii[2] ];
             _newTetras[2][1] = NodeGlobalID[ ii[7] ];
             _newTetras[2][2] = NodeGlobalID[ ii[8] ];
             _newTetras[2][3] = NodeGlobalID[ ii[6] ];

             _newTetras[3][0] = NodeGlobalID[ ii[1] ];
             _newTetras[3][1] = NodeGlobalID[ ii[8] ];
             _newTetras[3][2] = NodeGlobalID[ ii[3] ];
             _newTetras[3][3] = NodeGlobalID[ ii[4] ];

             _newTetras[4][0] = NodeGlobalID[ ii[1] ];
             _newTetras[4][1] = NodeGlobalID[ ii[8] ];
             _newTetras[4][2] = NodeGlobalID[ ii[2] ];
             _newTetras[4][3] = NodeGlobalID[ ii[3] ];

             _newTetras[5][0] = NodeGlobalID[ ii[2] ];
             _newTetras[5][1] = NodeGlobalID[ ii[8] ];
             _newTetras[5][2] = NodeGlobalID[ ii[7] ];
             _newTetras[5][3] = NodeGlobalID[ ii[3] ];

             _newTetrasSize = 6;
            break;

    case 2 : // 2 diagonals passing through vertex 7
             _newTetras[0][0] = NodeGlobalID[ ii[1] ];
             _newTetras[0][1] = NodeGlobalID[ ii[5] ];
             _newTetras[0][2] = NodeGlobalID[ ii[6] ];
             _newTetras[0][3] = NodeGlobalID[ ii[7] ];

             _newTetras[1][0] = NodeGlobalID[ ii[1] ];
             _newTetras[1][1] = NodeGlobalID[ ii[4] ];
             _newTetras[1][2] = NodeGlobalID[ ii[8] ];
             _newTetras[1][3] = NodeGlobalID[ ii[7] ];

             _newTetras[2][0] = NodeGlobalID[ ii[1] ];
             _newTetras[2][1] = NodeGlobalID[ ii[8] ];
             _newTetras[2][2] = NodeGlobalID[ ii[5] ];
             _newTetras[2][3] = NodeGlobalID[ ii[7] ];

             _newTetras[3][0] = NodeGlobalID[ ii[1] ];
             _newTetras[3][1] = NodeGlobalID[ ii[2] ];
             _newTetras[3][2] = NodeGlobalID[ ii[3] ];
             _newTetras[3][3] = NodeGlobalID[ ii[6] ];

             _newTetras[4][0] = NodeGlobalID[ ii[1] ];
             _newTetras[4][1] = NodeGlobalID[ ii[4] ];
             _newTetras[4][2] = NodeGlobalID[ ii[7] ];
             _newTetras[4][3] = NodeGlobalID[ ii[3] ];

             _newTetras[5][0] = NodeGlobalID[ ii[1] ];
             _newTetras[5][1] = NodeGlobalID[ ii[7] ];
             _newTetras[5][2] = NodeGlobalID[ ii[6] ];
             _newTetras[5][3] = NodeGlobalID[ ii[3] ];

             _newTetrasSize = 6;
            break;

    case 3 : // 3 diagonals passing through vertex 7
             _newTetras[0][0] = NodeGlobalID[ ii[1] ];
             _newTetras[0][1] = NodeGlobalID[ ii[3] ];
             _newTetras[0][2] = NodeGlobalID[ ii[4] ];
             _newTetras[0][3] = NodeGlobalID[ ii[7] ];

             _newTetras[1][0] = NodeGlobalID[ ii[1] ];
             _newTetras[1][1] = NodeGlobalID[ ii[4] ];
             _newTetras[1][2] = NodeGlobalID[ ii[8] ];
             _newTetras[1][3] = NodeGlobalID[ ii[7] ];

             _newTetras[2][0] = NodeGlobalID[ ii[1] ];
             _newTetras[2][1] = NodeGlobalID[ ii[8] ];
             _newTetras[2][2] = NodeGlobalID[ ii[5] ];
             _newTetras[2][3] = NodeGlobalID[ ii[7] ];

             _newTetras[3][0] = NodeGlobalID[ ii[1] ];
             _newTetras[3][1] = NodeGlobalID[ ii[6] ];
             _newTetras[3][2] = NodeGlobalID[ ii[7] ];
             _newTetras[3][3] = NodeGlobalID[ ii[5] ];

             _newTetras[4][0] = NodeGlobalID[ ii[2] ];
             _newTetras[4][1] = NodeGlobalID[ ii[6] ];
             _newTetras[4][2] = NodeGlobalID[ ii[7] ];
             _newTetras[4][3] = NodeGlobalID[ ii[1] ];

             _newTetras[5][0] = NodeGlobalID[ ii[2] ];
             _newTetras[5][1] = NodeGlobalID[ ii[7] ];
             _newTetras[5][2] = NodeGlobalID[ ii[3] ];
             _newTetras[5][3] = NodeGlobalID[ ii[1] ];

             _newTetrasSize = 6;
            break;
    default: CFout << "Problem in HEXAHEDRON Splitting" << "\n";
             for (CFuint i=0; i < 8 ; ++i) CFout << "Hexa: " << hexa[i] << "\n";
             cf_assert(0);
            break;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint CellSplitter3DFVM::nbNewTetras(CFuint iHexa)
{
  // Hexa in the article notation (that's why '+1')
  vector<CFuint> NodeGlobalID( nodesInHexa + 1 );
  // Indirections
  vector<CFuint> ii( nodesInHexa + 1 );

  // assign node numbering of the hexahedron
  for (CFuint i=0; i<nodesInHexa ; ++i){
    NodeGlobalID[i+1] = _oldElemNode(iHexa, i);
  }

  _SmallestFaceID = -1;
  getSmallestID(nodesInHexa, NodeGlobalID);
  cf_assert( _SmallestFaceID >= 0 );

  // create indirections = rotate hexahedron in order to have smallest index in lower-front-left position
  switch( _SmallestFaceID ) {
    case 1 : ii[1] = 1;
             ii[2] = 2;
             ii[3] = 3;
             ii[4] = 4;
             ii[5] = 5;
             ii[6] = 6;
             ii[7] = 7;
             ii[8] = 8;
             break;

    case 2 : ii[1] = 2;
             ii[2] = 1;
             ii[3] = 5;
             ii[4] = 6;
             ii[5] = 3;
             ii[6] = 4;
             ii[7] = 8;
             ii[8] = 7;
             break;

    case 3 : ii[1] = 3;
             ii[2] = 2;
             ii[3] = 6;
             ii[4] = 7;
             ii[5] = 4;
             ii[6] = 1;
             ii[7] = 5;
             ii[8] = 8;
             break;

    case 4 : ii[1] = 4;
             ii[2] = 1;
             ii[3] = 2;
             ii[4] = 3;
             ii[5] = 8;
             ii[6] = 5;
             ii[7] = 6;
             ii[8] = 7;
             break;

    case 5 : ii[1] = 5;
             ii[2] = 1;
             ii[3] = 4;
             ii[4] = 8;
             ii[5] = 6;
             ii[6] = 2;
             ii[7] = 3;
             ii[8] = 7;
             break;

    case 6 : ii[1] = 6;
             ii[2] = 2;
             ii[3] = 1;
             ii[4] = 5;
             ii[5] = 7;
             ii[6] = 3;
             ii[7] = 4;
             ii[8] = 8;
             break;

    case 7 : ii[1] = 7;
             ii[2] = 3;
             ii[3] = 2;
             ii[4] = 6;
             ii[5] = 8;
             ii[6] = 4;
             ii[7] = 1;
             ii[8] = 5;
             break;

    case 8 : ii[1] = 8;
             ii[2] = 4;
             ii[3] = 3;
             ii[4] = 7;
             ii[5] = 5;
             ii[6] = 1;
             ii[7] = 2;
             ii[8] = 6;
             break;

    default: CFerr << " When splitting HEXA into TETRA - face index out of range.\n";
             break;
  }

  // check whether there is any diagonal passing through corner (7)
  bool DiagThrough7 = false;

  if ( min( NodeGlobalID[ ii[2] ], NodeGlobalID[ ii[7] ] ) < min( NodeGlobalID[ ii[3] ], NodeGlobalID[ ii[6] ] ) )  DiagThrough7 = true;
  if ( min( NodeGlobalID[ ii[4] ], NodeGlobalID[ ii[7] ] ) < min( NodeGlobalID[ ii[3] ], NodeGlobalID[ ii[8] ] ) )  DiagThrough7 = true;
  if ( min( NodeGlobalID[ ii[5] ], NodeGlobalID[ ii[7] ] ) < min( NodeGlobalID[ ii[6] ], NodeGlobalID[ ii[8] ] ) )  DiagThrough7 = true;

  CFuint nbNewTetrasFromCurrentHexa;

  if( DiagThrough7 == false){  nbNewTetrasFromCurrentHexa = 5;  }
  else{                        nbNewTetrasFromCurrentHexa = 6;  }

  return nbNewTetrasFromCurrentHexa;
}
//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::getSmallestID(CFuint nbNodes, vector<CFuint> OneElement)
{
  // temporary variable for determining position of the smallest index
  vector<CFuint> Order( nbNodes + 1 );
  // local copy of node numbering (otherwise messing up numbering in global context)
  vector<CFuint> ID( nbNodes + 1 );

  for (CFuint i=0; i<nbNodes ; ++i){
    ID[i] = OneElement[i+1];
    Order[i] = i+1; // (plus 1) because of article numbering
  }

  for( CFuint i=nbNodes-1; i>0; i-- ){
    for( CFuint j=0; j<=i-1; j++){
      if( ID[j] > ID[j+1] ){
        CFuint tmpID = ID[j];
        ID[j]        = ID[j+1];
        ID[j+1]      = tmpID;

        CFuint tmpOrder = Order[j];
        Order[j]        = Order[j+1];
        Order[j+1]      = tmpOrder;
      }
    }
  }

  _SmallestFaceID = Order[0];
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::updateTRSData()
{
  CFAUTOTRACE;

  CFmeshReaderWriterSource& data = *(_data.get());

  SafePtr< vector<CFuint> > nbTRs = data.getNbTRs();
  SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR = data.getNbGeomEntsPerTR();

  SafePtr< Table<CFuint> > elementNode  = data.getElementNode();
  const CFuint nbElem = data.getNbElements();

  // loop over boundaries (TRS)
  for(CFuint iTRS = 0; iTRS < data.getNbTRSs(); ++iTRS ) {
    // read connectivity of iTRS's boundary
    TRGeoConn& geoConn = data.getTRGeoConn(iTRS);
    cf_assert((*nbTRs)[iTRS] == geoConn.size());

    // loop over all TR inside current TRS
    for(CFuint iTR = 0; iTR < (*nbTRs)[iTRS]; ++iTR ) {

      // duplicate old connectivity list
            GeoConn& geoC    = geoConn[iTR];
      const GeoConn& oldGeoC = geoConn[iTR];
      const CFuint oldNbGeos = geoC.size();
      GeoConn oldConn = oldGeoC;
      // computes how many new elements will be created
      CFuint newNbGeos = 0;
      for( CFuint iGeo=0; iGeo < oldNbGeos; ++iGeo ){
        std::valarray<CFuint>& oldNodeCon   = (oldConn[iGeo].first);

        if( oldNodeCon.size() == 3 ) // triangle
          newNbGeos ++;
        if( oldNodeCon.size() == 4 ) // quad
          newNbGeos += 2;
      }
      // resize (and clean) original list geoC
      geoC.resize(newNbGeos);

///THIS IS WRONG FOR MULTIELEMENT GRIDS
//       // fill first part with old data
//       for( CFuint oldID=0; oldID < oldNbGeos; ++oldID ){
//         geoC[ oldID ] = oldGeoC[ oldID ];
//       }

      // loop directly over elements
      CFuint shift = 0;
      for (CFuint iGeo = 0; iGeo < oldNbGeos; ++iGeo) {
        std::valarray<CFuint>& oldNodeCon   = (oldConn[iGeo].first);

        if ( oldNodeCon.size()==3 ) {
          // Triangle - connectivity must be just copied
          geoC[ iGeo ] = oldGeoC[ iGeo ];
        }

        if ( oldNodeCon.size()==4 ) {
          // Quad - will be splitted into 2 triangles (along smallestID node)
          vector<CFuint> tempQuad;
          tempQuad.resize(nodesInQuad);
          for(CFuint quadID = 0; quadID < nodesInQuad; ++quadID){
            tempQuad[quadID] = oldNodeCon[quadID];
          }
          splitQuads(tempQuad);

          // Set the new node connectivity
          std::valarray<CFuint>& nodeCon1 = geoC[ iGeo ].first;
          std::valarray<CFuint>& nodeCon2 = geoC[ oldNbGeos + shift ].first;
          nodeCon1.resize(3);
          nodeCon2.resize(3);

          for(CFuint nodeID = 0; nodeID < 3; ++nodeID) {
            nodeCon1[nodeID] = _newTriags[0][nodeID];
            nodeCon2[nodeID] = _newTriags[1][nodeID];
          }

          // Set the new state connectivity
          std::valarray<CFuint>& stateCon1 = geoC[ iGeo ].second;
          std::valarray<CFuint>& stateCon2 = geoC[ oldNbGeos + shift ].second;
          stateCon1.resize(1);
          stateCon2.resize(1);

          // localize the state ID - must loop over all elements and compare nodes
          for(CFuint iNewTriang=0; iNewTriang < 2; iNewTriang++){

            const CFuint N0 = _newTriags[ iNewTriang ][ 0 ];
            const CFuint N1 = _newTriags[ iNewTriang ][ 1 ];
            const CFuint N2 = _newTriags[ iNewTriang ][ 2 ];

            for(CFuint iElem=0; iElem < nbElem; iElem++){ // loop over all elements
              for(CFuint iN0=0; iN0 < 4; iN0++){ // loop over all nodes of element
                if( N0 == (*elementNode)(iElem,iN0) ){ // agreement on first node

                  for(CFuint iN1=0; iN1 < 4; iN1++){
                    if( N1 == (*elementNode)(iElem,iN1) ){ // agreement on second node

                      for(CFuint iN2=0; iN2 < 4; iN2++){
                        if( N2 == (*elementNode)(iElem,iN2) ){ // agreement on third node - search completed

                          if( iNewTriang == 0 ) stateCon1[0] = iElem;
                          if( iNewTriang == 1 ) stateCon2[0] = iElem;
                      } }
                  } }
              } }
            }
          }
          shift++;
        }
        if( (oldNodeCon.size() !=3 ) && (oldNodeCon.size() !=4 ) ) {
          std::string msg = std::string("Can't split element on TRS, because it is neither Triangle nor Quad.");
          throw BadValueException (FromHere(),msg);
        }
      }

      (*nbGeomEntsPerTR)[iTRS][iTR] = geoC.size();

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellSplitter3DFVM::splitQuads(vector<CFuint> quad)
{
  const CFuint NodeA = quad[0];
  const CFuint NodeB = quad[1];
  const CFuint NodeC = quad[2];
  const CFuint NodeD = quad[3];

  // smallest ID at node A
  if( (NodeA<NodeB) && (NodeA<NodeC) && (NodeA<NodeD) ){
    _newTriags[0][0] = NodeA;
    _newTriags[0][1] = NodeB;
    _newTriags[0][2] = NodeC;

    _newTriags[1][0] = NodeC;
    _newTriags[1][1] = NodeD;
    _newTriags[1][2] = NodeA;
  }

  // smallest ID at node B
  if( (NodeB<NodeA) && (NodeB<NodeC) && (NodeB<NodeD) ){
    _newTriags[0][0] = NodeB;
    _newTriags[0][1] = NodeC;
    _newTriags[0][2] = NodeD;

    _newTriags[1][0] = NodeD;
    _newTriags[1][1] = NodeA;
    _newTriags[1][2] = NodeB;
  }

  // smallest ID at node C
  if( (NodeC<NodeA) && (NodeC<NodeB) && (NodeC<NodeD) ){
    _newTriags[0][0] = NodeC;
    _newTriags[0][1] = NodeD;
    _newTriags[0][2] = NodeA;

    _newTriags[1][0] = NodeA;
    _newTriags[1][1] = NodeB;
    _newTriags[1][2] = NodeC;
  }

  // smallest ID at node D
  if( (NodeD<NodeA) && (NodeD<NodeB) && (NodeD<NodeC) ){
    _newTriags[0][0] = NodeD;
    _newTriags[0][1] = NodeA;
    _newTriags[0][2] = NodeB;

    _newTriags[1][0] = NodeB;
    _newTriags[1][1] = NodeC;
    _newTriags[1][2] = NodeD;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace CFmeshCellSplitter

  } // end of namespace IO

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
