// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>

#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/Stopwatch.hh"
#include "Common/CFMap.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/Table.hh"
#include "Gambit2CFmesh/Gambit2CFmesh.hh"
#include "Gambit2CFmesh/Gambit2CFmeshConverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Gambit2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Gambit2CFmeshConverter,
               MeshFormatConverter,
               Gambit2CFmeshModule,
               1>
gambit2CFmeshConverterProvider("Gambit2CFmesh");

//////////////////////////////////////////////////////////////////////////////

Gambit2CFmeshConverter::Gambit2CFmeshConverter (const std::string& name)
: MeshFormatConverter(name),
  _dimension(0),
  _nbUpdatableNodes(0),
  _nbCells(0),
  _nbFaces(0),
  _nbPatches(0),
  _nbSuperPatches(0),
  _elementType(0),
  _patch(0),
  _superPatch(0),
  _update(static_cast<CFuint>(0),static_cast<CFuint>(0)),
  _coordinate(),
  _variables(CFNULL),
  _isWithSolution(false),
  _isFileRead(false),
  _nbUpdatableStates(_nbUpdatableNodes),
  _nodesPerElemTypeTable(7),
  _inFieldSP(0),
  _nbGroups(0),
  _nbVelComponents(0),
  _isHighOrder(false),
  _orderChangeMap(CFNULL),
  _faceDefinionTable(0),
  _tablePattern(0),
  _gambitElemTypeID(0),
  _mapFaces(),
  _mapNodeOrderHexa(),
  _mapNodeOrderPyram()
{
  /// Build the nbNodes per ElemTypeTable
  _nodesPerElemTypeTable[0] = 2; // Edge
  _nodesPerElemTypeTable[1] = 3; // Triag
  _nodesPerElemTypeTable[2] = 4; // Quad
  _nodesPerElemTypeTable[3] = 4; // Tetra
  _nodesPerElemTypeTable[4] = 8; // Hexa
  _nodesPerElemTypeTable[5] = 6; // Prism
  _nodesPerElemTypeTable[6] = 5; // Pyram

  /// Tables of edge or face definitions

  /**
  * Edge - 1x2 table
  */
  _faceDefinionTable.push_back(new Table<CFuint>(1,2));

  (*_faceDefinionTable[0])(0,0) = 0;
  (*_faceDefinionTable[0])(0,1) = 1;

  /**
  * Triag - 3x2 table
  */
  _faceDefinionTable.push_back(new Table<CFuint>(3,2));

  (*_faceDefinionTable[1])(0,0) = 0;
  (*_faceDefinionTable[1])(0,1) = 1;

  (*_faceDefinionTable[1])(1,0) = 1;
  (*_faceDefinionTable[1])(1,1) = 2;

  (*_faceDefinionTable[1])(2,0) = 2;
  (*_faceDefinionTable[1])(2,1) = 0;

  /**
  * Quag - 4x2 table
  */
  _faceDefinionTable.push_back(new Table<CFuint>(4,2));

  (*_faceDefinionTable[2])(0,0) = 0;
  (*_faceDefinionTable[2])(0,1) = 1;

  (*_faceDefinionTable[2])(1,0) = 1;
  (*_faceDefinionTable[2])(1,1) = 2;

  (*_faceDefinionTable[2])(2,0) = 2;
  (*_faceDefinionTable[2])(2,1) = 3;

  (*_faceDefinionTable[2])(3,0) = 3;
  (*_faceDefinionTable[2])(3,1) = 0;

  /**
  * Tetra - 4x3 table
  */
  _faceDefinionTable.push_back(new Table<CFuint>(4,3));

  (*_faceDefinionTable[3])(0,0) = 0;
  (*_faceDefinionTable[3])(0,1) = 2;
  (*_faceDefinionTable[3])(0,2) = 1;

  (*_faceDefinionTable[3])(1,0) = 0;
  (*_faceDefinionTable[3])(1,1) = 1;
  (*_faceDefinionTable[3])(1,2) = 3;

  (*_faceDefinionTable[3])(2,0) = 1;
  (*_faceDefinionTable[3])(2,1) = 2;
  (*_faceDefinionTable[3])(2,2) = 3;

  (*_faceDefinionTable[3])(3,0) = 0;
  (*_faceDefinionTable[3])(3,1) = 3;
  (*_faceDefinionTable[3])(3,2) = 2;

  /**
  * Hexa - 6x4 table
  */
  _faceDefinionTable.push_back(new Table<CFuint>(6,4));

  (*_faceDefinionTable[4])(0,0) = 0;
  (*_faceDefinionTable[4])(0,1) = 3;
  (*_faceDefinionTable[4])(0,2) = 2;
  (*_faceDefinionTable[4])(0,3) = 1;

  (*_faceDefinionTable[4])(1,0) = 1;
  (*_faceDefinionTable[4])(1,1) = 2;
  (*_faceDefinionTable[4])(1,2) = 6;
  (*_faceDefinionTable[4])(1,3) = 5;

  (*_faceDefinionTable[4])(2,0) = 4;
  (*_faceDefinionTable[4])(2,1) = 5;
  (*_faceDefinionTable[4])(2,2) = 6;
  (*_faceDefinionTable[4])(2,3) = 7;

  (*_faceDefinionTable[4])(3,0) = 3;
  (*_faceDefinionTable[4])(3,1) = 0;
  (*_faceDefinionTable[4])(3,2) = 4;
  (*_faceDefinionTable[4])(3,3) = 7;

  (*_faceDefinionTable[4])(4,0) = 2;
  (*_faceDefinionTable[4])(4,1) = 3;
  (*_faceDefinionTable[4])(4,2) = 7;
  (*_faceDefinionTable[4])(4,3) = 6;

  (*_faceDefinionTable[4])(5,0) = 0;
  (*_faceDefinionTable[4])(5,1) = 1;
  (*_faceDefinionTable[4])(5,2) = 5;
  (*_faceDefinionTable[4])(5,3) = 4;

  /**
  * Prism
  */

  /// Definition of table structure for Prism
  _tablePattern.push_back(std::valarray<CFuint>(5));
  _tablePattern[0][0] = 4;
  _tablePattern[0][1] = 4;
  _tablePattern[0][2] = 4;
  _tablePattern[0][3] = 3;
  _tablePattern[0][4] = 3;

  /// Initialization of the table for Prism
  _faceDefinionTable.push_back(new Table<CFuint>(_tablePattern[0]));

  (*_faceDefinionTable[5])(0,0) = 0;
  (*_faceDefinionTable[5])(0,1) = 1;
  (*_faceDefinionTable[5])(0,2) = 4;
  (*_faceDefinionTable[5])(0,3) = 3;

  (*_faceDefinionTable[5])(1,0) = 1;
  (*_faceDefinionTable[5])(1,1) = 2;
  (*_faceDefinionTable[5])(1,2) = 5;
  (*_faceDefinionTable[5])(1,3) = 4;

  (*_faceDefinionTable[5])(2,0) = 2;
  (*_faceDefinionTable[5])(2,1) = 0;
  (*_faceDefinionTable[5])(2,2) = 3;
  (*_faceDefinionTable[5])(2,3) = 5;

  (*_faceDefinionTable[5])(3,0) = 0;
  (*_faceDefinionTable[5])(3,1) = 2;
  (*_faceDefinionTable[5])(3,2) = 1;

  (*_faceDefinionTable[5])(4,0) = 3;
  (*_faceDefinionTable[5])(4,1) = 4;
  (*_faceDefinionTable[5])(4,2) = 5;

  /**
  * Pyram
  */

  /// Definition of table structure for Pyram
  _tablePattern.push_back(std::valarray<CFuint>(5));
  _tablePattern[1][0] = 4;
  _tablePattern[1][1] = 3;
  _tablePattern[1][2] = 3;
  _tablePattern[1][3] = 3;
  _tablePattern[1][4] = 3;

  /// Initialization of the table for Pyram
  _faceDefinionTable.push_back(new Table<CFuint>(_tablePattern[1]));

  (*_faceDefinionTable[6])(0,0) = 0;
  (*_faceDefinionTable[6])(0,1) = 3;
  (*_faceDefinionTable[6])(0,2) = 2;
  (*_faceDefinionTable[6])(0,3) = 1;

  (*_faceDefinionTable[6])(1,0) = 0;
  (*_faceDefinionTable[6])(1,1) = 1;
  (*_faceDefinionTable[6])(1,2) = 4;

  (*_faceDefinionTable[6])(2,0) = 1;
  (*_faceDefinionTable[6])(2,1) = 2;
  (*_faceDefinionTable[6])(2,2) = 4;

  (*_faceDefinionTable[6])(3,0) = 2;
  (*_faceDefinionTable[6])(3,1) = 3;
  (*_faceDefinionTable[6])(3,2) = 4;

  (*_faceDefinionTable[6])(4,0) = 3;
  (*_faceDefinionTable[6])(4,1) = 0;
  (*_faceDefinionTable[6])(4,2) = 4;

  /**
  * Initialization of _gambitElemTypeID vector
  */

  _gambitElemTypeID.push_back(make_pair(1,2)); // Edge
  _gambitElemTypeID.push_back(make_pair(3,3)); // Triag
  _gambitElemTypeID.push_back(make_pair(2,4)); // Quad
  _gambitElemTypeID.push_back(make_pair(6,4)); // Tetra
  _gambitElemTypeID.push_back(make_pair(4,8)); // Hexa
  _gambitElemTypeID.push_back(make_pair(5,6)); // Prism
  _gambitElemTypeID.push_back(make_pair(7,5)); // Pyram

  /**
  * _gambitElemTypeID to _faceDefinionTable map
  */
  for (CFuint i = 0; i < _nodesPerElemTypeTable.size(); ++i) {
    _mapFaces[_gambitElemTypeID[i]] = _faceDefinionTable[i];
  }

  /**
  * map of elemenet node orders Gambit -> CF
  */

  /// Hexa
  _mapNodeOrderHexa[0] = 3;
  _mapNodeOrderHexa[1] = 2;
  _mapNodeOrderHexa[2] = 7;
  _mapNodeOrderHexa[3] = 6;
  _mapNodeOrderHexa[4] = 0;
  _mapNodeOrderHexa[5] = 1;
  _mapNodeOrderHexa[6] = 4;
  _mapNodeOrderHexa[7] = 5;

  /// Pyram
  _mapNodeOrderPyram[0] = 0;
  _mapNodeOrderPyram[1] = 1;
  _mapNodeOrderPyram[2] = 3;
  _mapNodeOrderPyram[3] = 2;
  _mapNodeOrderPyram[4] = 4;
}

////////////////////////////////////////////////////////////////////////////// 0

Gambit2CFmeshConverter::~Gambit2CFmeshConverter()
{
  deletePtr(_variables);
  deletePtr(_orderChangeMap);
  for (CFuint i = 0; i < _faceDefinionTable.size(); ++i) deletePtr(_faceDefinionTable[i]);
}

////////////////////////////////////////////////////////////////////////////// 1

void Gambit2CFmeshConverter::checkFormat(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
  path meshFile = change_extension(filepath, getOriginExtension());

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(meshFile);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;

  // Check for GAMBIT
  // The second line must be with 4 words = "** GAMBIT NEUTRAL FILE"
  getGambitWordsFromLine(fin,line,lineNb,words);  // L.1
  getGambitWordsFromLine(fin,line,lineNb,words);  // L.2

  if(words.size() != 4) callGambitFileError("Malformed file format - not a GAMBIT NEUTRAL FILE format",lineNb,meshFile);
  if((words[0] != "**")&&(words[1] != "GAMBIT")&&(words[2] != "NEUTRAL")&&(words[3] != "FILE")) callGambitFileError("Malformed file format - not a GAMBIT NEUTRAL FILE format ",lineNb,meshFile);

  // Skip lines
  getGambitWordsFromLine(fin,line,lineNb,words);  // L.3
  getGambitWordsFromLine(fin,line,lineNb,words);  // L.4
  getGambitWordsFromLine(fin,line,lineNb,words);  // L.5

  // Check nb Elements,
  // nb nodes, nb Boundary Faces
  // nb patches
  getGambitWordsFromLine(fin,line,lineNb,words);  // L.6
  if(words.size() > 6) callGambitFileError("Wrong number of parameters.",lineNb,meshFile);
  if ((words[0] != "NUMNP")&&(words[1] != "NELEM")&&(words[2] != "NGRPS")&&(words[3] != "NBSETS")&&(words[4] != "NDFCD")&&(words[5] != "NDFVL")) callGambitFileError("Malformed file format - not GAMBIT NEUTRAL FILE format ",lineNb,meshFile);

  getGambitWordsFromLine(fin,line,lineNb,words);  // L.7
  if(words.size() != 6) callGambitFileError("Wrong number of parameters.",lineNb,meshFile);

  CFint const nbNodes = StringOps::from_str<CFint>(words[0]);
  if(nbNodes <= 0) callGambitFileError("Negative number of Nodes.",lineNb,meshFile);

  CFint const nbCells = StringOps::from_str<CFint>(words[1]);
  if(nbCells < 0) callGambitFileError("Negative number of Elements.",lineNb,meshFile);

  CFint const nbGroups = StringOps::from_str<CFint>(words[2]);
  if(nbGroups < 0) callGambitFileError("Negative number of Groups.",lineNb,meshFile);

  CFint const nbSuperPatches = StringOps::from_str<CFint>(words[3]);
  if(nbSuperPatches < 0) callGambitFileError("Negative number of Boundary Condition Sets.",lineNb,meshFile);

  CFint const dimension = StringOps::from_str<CFint>(words[4]);
  if((dimension < 2) || (dimension > 3)) callGambitFileError("Wrong number of Coordinate Directions.",lineNb,meshFile);

  CFint const nbVelComponents = StringOps::from_str<CFint>(words[5]);
  if((nbVelComponents < 2) || (nbVelComponents > 3)) callGambitFileError("Wrong number of Velocity Components.",lineNb,meshFile);

  // Close file
  fhandle->close();
}

////////////////////////////////////////////////////////////////////////////// 2

void Gambit2CFmeshConverter::callGambitFileError(const std::string& msg,
               const CFuint& line,
               const boost::filesystem::path& filepath)
{
  std::string fullMessage("Line ");
  fullMessage += StringOps::to_str(line);
  fullMessage += " File ";
  fullMessage += filepath.string();
  fullMessage += ". ";
  fullMessage += msg;
  throw BadFormatException (FromHere(),fullMessage);
}

////////////////////////////////////////////////////////////////////////////// 3

void Gambit2CFmeshConverter::getGambitWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

////////////////////////////////////////////////////////////////////////////// 3

void Gambit2CFmeshConverter::getGroupWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;

  vector<std::string> wordsLocal;
  const CFuint dataInLine = line.length() / 8; // 8 characters/data
  for(CFuint i=0; i<dataInLine; i++)
  {
    wordsLocal.push_back(line.substr(i*8,8));
  }
  words = wordsLocal;
}

////////////////////////////////////////////////////////////////////////////// 3

void Gambit2CFmeshConverter::getFirstConnectivityWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;

  vector<std::string> wordsLocal;
  wordsLocal.push_back(line.substr(0,8));  // element/cell ID
  wordsLocal.push_back(line.substr(9,2));  // element/cell type
  wordsLocal.push_back(line.substr(12,2)); // number of nodes in element/cell
  const CFuint dataInLine = (line.length()-15) / 8; // 8 characters/data after position 15
  for(CFuint i=0; i<dataInLine; i++) // read node IDs
  {
    wordsLocal.push_back(line.substr(15+i*8,8));
  }
  words = wordsLocal;
}

////////////////////////////////////////////////////////////////////////////// 3

void Gambit2CFmeshConverter::getRestConnectivityWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;

  vector<std::string> wordsLocal;
  const CFuint dataInLine = (line.length()-15) / 8; // 8 characters/data after position 15
  for(CFuint i=0; i<dataInLine; i++) // read node IDs
  {
    wordsLocal.push_back(line.substr(15+i*8,8));
  }
  words = wordsLocal;
}

////////////////////////////////////////////////////////////////////////////// 4

void Gambit2CFmeshConverter::readParameters(ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb)
{
  CFAUTOTRACE;

  std::string line;
  vector<std::string> words;

  // Jump to 7-th line of Gambit neutral file
  for(CFuint i=1; i<=6; i++)  getGambitWordsFromLine(fin,line,lineNb,words); // L.1 - L.6

  // Reads line 7.
  getGambitWordsFromLine(fin,line,lineNb,words); // L.7

  // Stores data from line 7
  _nbUpdatableNodes = StringOps::from_str<CFuint>(words[0]);
  _nbCells          = StringOps::from_str<CFuint>(words[1]);
  _nbGroups         = StringOps::from_str<CFuint>(words[2]);
  _nbSuperPatches   = StringOps::from_str<CFuint>(words[3]);
  _dimension        = StringOps::from_str<CFuint>(words[4]);
  _nbVelComponents  = StringOps::from_str<CFuint>(words[5]);

  _nbPatches = _nbSuperPatches;
}

////////////////////////////////////////////////////////////////////////////// 5

void Gambit2CFmeshConverter::readNodes(ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb)
{
  CFAUTOTRACE;

  std::string line;
  vector<std::string> words;

  // Table of node coordinaates
  _coordinate.resize(_nbUpdatableNodes*_dimension);

  // Reads coordinates from file and stores them in table
  CFuint prevNodeID = 0;
  for (CFuint i = 0; i < _nbUpdatableNodes; ++i){
    getGambitWordsFromLine(fin,line,lineNb,words);

    // First node ID should be 1 and the IDs should increase by 1 (in this mesh converter)
    CFuint currNodeID = StringOps::from_str<CFuint>(words[0]);
    if( currNodeID !=  prevNodeID + 1) callGambitFileError("Wrong ordering of Nodes. Node IDs must follow each other.",lineNb,filepath);
    prevNodeID = currNodeID;

    for (CFuint j = 0; j < _dimension; ++j){
      _coordinate[i*_dimension+j] = StringOps::from_str<CFreal>(words[j+1]);
    }
  }

  // Checks if the ENDSECTION is reached
  getGambitWordsFromLine(fin,line,lineNb,words);
  if(words[0] != "ENDOFSECTION") callGambitFileError("Wrong number of Nodes. ENDOFSECTION not reached.",lineNb,filepath);
}

////////////////////////////////////////////////////////////////////////////// 5

void Gambit2CFmeshConverter::readGroups(ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb)
{
  CFAUTOTRACE;

  std::string line;
  vector<std::string> words;

  _groupCellIDs.resize(_nbGroups);
  _groupNames.resize(_nbGroups);
  for(CFuint iGroup=0; iGroup < _nbGroups; ++iGroup){
    getGambitWordsFromLine(fin,line,lineNb,words);
    if((words[0] != "ELEMENT") || (words[1] != "GROUP")) callGambitFileError("Wrong order of section, no group found.",lineNb,filepath);
    getGambitWordsFromLine(fin,line,lineNb,words);
    if(words[0] != "GROUP:") callGambitFileError("Wrong order of section, group not found.",lineNb,filepath);

    const CFuint groupID = StringOps::from_str<CFuint>(words[1]);
    if(groupID != (iGroup+1)) callGambitFileError("Wrong order of section, group description not found.",lineNb,filepath);

    cf_assert(words[2] == "ELEMENTS:");
    const CFuint nbGroupElements = StringOps::from_str<CFuint>(words[3]);
    _groupCellIDs[iGroup].resize(nbGroupElements);

    getGambitWordsFromLine(fin,line,lineNb,words);
    if(words.size() != 1) callGambitFileError("Malformed file. Group Name could not be found.",lineNb,filepath);
    _groupNames[iGroup] = words[0];

    //read the solver dependent flag (usually 0)
    getGambitWordsFromLine(fin,line,lineNb,words);
    cf_assert(words.size() == 1);

    CFuint iElem = 0;
    while(iElem < nbGroupElements){
      getGroupWordsFromLine(fin,line,lineNb,words);
      const CFuint nbElementsInLine = words.size();
      // Reads ElementIDs from file and stores them
      for (CFuint i = 0; i < nbElementsInLine; ++i){
        _groupCellIDs[iGroup][iElem] = StringOps::from_str<CFuint>(words[i]);
        ++iElem;
      }
    }
//    cf_assert(iElem == nbGroupElements);

    // Checks if the ENDSECTION is reached
    getGambitWordsFromLine(fin,line,lineNb,words);
    if(words[0] != "ENDOFSECTION") callGambitFileError("Wrong number of Elements in group. ENDOFSECTION not reached.",lineNb,filepath);

  }


}


////////////////////////////////////////////////////////////////////////////// 6

void Gambit2CFmeshConverter::readConnectivity(ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb)
{
  CFAUTOTRACE;

  // Variable declaration
  std::string line;
  vector<std::string> words;

  CFLogDebugMin( "Gambit2CFmeshConverter::readConnectivity() => start\n" );

  CFuint maxNbElementTypes = _nodesPerElemTypeTable.size(); // Max possible number of considered elements
  CFuint nbElementTypes = 0; // Number of elementTypes with non-zero elements
  CFuint gambitElementTypeID = 0; // Number of element type in Gambit neutral format
  CFuint nbNodesPerElemType = 0; // Number of nodes per element type
  CFuint maxNbNodesPerCell = 0; // Maximum number of nodes per cell
  CFuint nbWordsOnLine = 0; // Number of words read on a line
  CFuint idx = 0;

  std::vector<CFuint> nbElemPerType(maxNbElementTypes); // Number of elements per type
  for (CFuint i=0; i < maxNbElementTypes; ++i) {nbElemPerType[i] = 0;} // Initialization of nbElemPerType vector

  std::vector<CFuint> vectorTypeID(_nbCells); // Vector of typeID
  for (CFuint i=0; i < _nbCells; ++i) {vectorTypeID[i] = 0;} // Initialization of vectorTypeID vector

  // Temporary table connectivity
  if (_isHighOrder) maxNbNodesPerCell = 27;
  else maxNbNodesPerCell = 8;

  std::valarray<CFuint> columnPattern(_nbCells);
  columnPattern = maxNbNodesPerCell;
  _tableConnectivity.resize(columnPattern);

  // Reads elements
  CFuint prevCellID = 0;
  for (CFuint i = 0; i < _nbCells; ++i)
  {
   // Reads first element data
   getFirstConnectivityWordsFromLine(fin,line,lineNb,words);
   nbWordsOnLine = words.size();

   // First cell/element ID should be 1 and the IDs should increase by 1 (in this mesh converter)
   CFuint currCellID = StringOps::from_str<CFuint>(words[0]);
   if( currCellID !=  prevCellID + 1) callGambitFileError("Wrong ordering of Cells. Cell IDs must follow each other.",lineNb,filepath);
   prevCellID = currCellID;

   // Reads number of element type in Gambit neutral format
   gambitElementTypeID = StringOps::from_str<CFuint>(words[1]);

   // Reads number of nodes per element type
   nbNodesPerElemType = StringOps::from_str<CFuint>(words[2]);

   vector<CFuint> nodesOfElement(nbNodesPerElemType);
   for(CFuint iNode = 0; iNode < (nbWordsOnLine-3); ++iNode){
     nodesOfElement[iNode] = (StringOps::from_str<CFuint>(words[iNode+3]))-1;
   }

   CFuint lastIdx = nbWordsOnLine-3;
   bool isAchieved = false;
   if(nbWordsOnLine != nbNodesPerElemType+3)
   {
    isAchieved = false;
    while(isAchieved != true){
         getRestConnectivityWordsFromLine(fin,line,lineNb,words);
         for(CFuint iNode = 0; iNode < words.size(); ++iNode){
           nodesOfElement[iNode+lastIdx] = (StringOps::from_str<CFuint>(words[iNode]))-1;
         }
         lastIdx += words.size();
         if(lastIdx == nbNodesPerElemType) isAchieved = true;
    }
   }

   switch(gambitElementTypeID)
   {
    // Gambit numbering
    case 1: // Edge
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 2: vectorTypeID[i] = 0; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
                     for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
                       _tableConnectivity(i,j) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Edge type.",lineNb,filepath);
      }
      break;
    case 2: // Quad
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 4: vectorTypeID[i] = 2; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
               for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
                       _tableConnectivity(i,j) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Quad type.",lineNb,filepath);
      }
      break;
    case 3: // Triag
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 3: vectorTypeID[i] = 1; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
               for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
                       _tableConnectivity(i,j) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Triag type.",lineNb,filepath);
      }
      break;
    case 4: // Hexa
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 8: vectorTypeID[i] = 4; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
               for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
           idx = _mapNodeOrderHexa.find(j)->second;
                       _tableConnectivity(i,idx) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Hexa type.",lineNb,filepath);
      }
      break;
    case 5: // Prism
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 6: vectorTypeID[i] = 5; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
               for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
                       _tableConnectivity(i,j) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Prism type.",lineNb,filepath);
      }
      break;
    case 6: // Tetra
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 4: vectorTypeID[i] = 3; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
               for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
                       _tableConnectivity(i,j) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Tetra type.",lineNb,filepath);
      }
      break;
    case 7: // Pyram
            switch(nbNodesPerElemType)
      {
       // Number of nodes
       case 5: vectorTypeID[i] = 6; // Stores CF element type number in vector
               ++nbElemPerType[vectorTypeID[i]];
               for (CFuint j = 0; j < nbNodesPerElemType; ++j) {
           idx = _mapNodeOrderPyram.find(j)->second;
                       _tableConnectivity(i,idx) = nodesOfElement[j];
         }
         break;
       default: // Other node number
                callGambitFileError("Number of nodes not implemented for Pyram type.",lineNb,filepath);
      }
      break;
    default: callGambitFileError("Not implemented Gambit element type.",lineNb,filepath);
   }
  }

  // Count the number of elementTypes with non-zero elements
  for (CFuint i=0; i < maxNbElementTypes; ++i){
    if (nbElemPerType[i] > 0) ++nbElementTypes;
  }

  // For each elementType present in the mesh, set the typeID
  // Read and allocate element type information
  _elementType.resize(nbElementTypes);

  // Sets element type ID with respect to CF element numbering
  CFuint *elemTypeOrder= new CFuint[maxNbElementTypes]; // Stores order of each element type in array of element types existing in mesh file
  for (CFuint i = 0; i < maxNbElementTypes; ++i) elemTypeOrder[i] = numeric_limits<CFuint>::max(); // Initialization of elemTypeOrder

  CFuint k=0;
  for (CFuint i = 0; i < maxNbElementTypes; ++i){
    if (nbElemPerType[i] > 0){
      _elementType[k].setTypeID(i); // Sets element type ID
      elemTypeOrder[i] = k; // Sets order of element type
      ++k;
    }
  }

  // For each elementType present in the mesh, set the nbNodesPerCell and NbCellsPerType
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    _elementType[k].setNbNodesPerCell(_nodesPerElemTypeTable[_elementType[k].getTypeID()]);
    _elementType[k].setNbCellsPerType(nbElemPerType[_elementType[k].getTypeID()]);
  }

  // Create storage space for the connectivity tables
  for (CFuint k = 0; k < nbElementTypes; ++k) {
    _elementType[k].createCellNodeConnectivity(&_tableConnectivity);
  }

  // Creates map connecting elements order to element type order
  createOrderChangeMap();

  CFuint *tableConnectivityLine= new CFuint[nbElementTypes]; // Counter of lines for each element type connectivity table
  for (CFuint i = 0; i < nbElementTypes; ++i) tableConnectivityLine[i] = 0; // Initialization of tableConnectivityLine

  // Copy tableConnectivity to _elementConnectivity
  for (CFuint i = 0; i < _nbCells; ++i) {
    nbNodesPerElemType = _nodesPerElemTypeTable[vectorTypeID[i]];
    k = elemTypeOrder[vectorTypeID[i]];
    _elementType[k].setElementGlobalID(i,tableConnectivityLine[k]);
    getOrderChangeMap()(i,0) = k;
    getOrderChangeMap()(i,1) = tableConnectivityLine[k];

    ++tableConnectivityLine[k];
  }

  // Checks if the ENDSECTION is reached
  getGambitWordsFromLine(fin,line,lineNb,words);
  if(words[0] != "ENDOFSECTION") callGambitFileError("Wrong number of Elements/Cells. ENDOFSECTION not reached.",lineNb,filepath);

  // freeing memory
	delete [] elemTypeOrder;
  delete [] tableConnectivityLine;

  CFLogDebugMin(  "Gambit2CFmeshConverter::readConnectivity() => end\n" );
}

////////////////////////////////////////////////////////////////////////////// 7

void Gambit2CFmeshConverter::readBC(ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb, const CFuint& patchCode)
{
  CFAUTOTRACE;

  CFLogDebugMin(  "Gambit2CFmeshConverter::readBC() => start\n" );

  // Variable declaration
  std::string line;
  vector<std::string> words;

  SPdata superPatch; // Temporary - in this method
  PatchGambit patch; // Temporary - in this method

  CFuint  firstFaceNode = 0;
  CFuint  secondFaceNode = 0;

  CFuint  gambitFaceID = 0;
  CFuint  gambitElemType = 0;
  CFuint  cellID = 0;
  CFuint  elemTypeOrder = 0;
  CFuint  lineInElemTypeTable = 0;
  CFuint  nbNodesPerElemType = 0;
  CFuint  nbFaceNodes = 0;
  CFuint  gambitFileNodeNumber = 0;
  CFuint  nodeOrderInElement = 0;
  CFuint  globalElemID = 0;
  Table<CFuint>*  faceDefinionTable = CFNULL;
  pair<CFuint, CFuint> gambitElemTypeID(0,0);

  CFuint bcType = 0; // Boundary condition type
  CFuint nbDataRecords = 0; // Number of Nodes or Elements/Cells in boundary condition

  // Reads header of BOUNDARY CONDITIONS section
  getGambitWordsFromLine(fin,line,lineNb,words);

  bcType = StringOps::from_str<CFuint>(words[1]);
  nbDataRecords = StringOps::from_str<CFuint>(words[2]);

  superPatch.setSuperPatchName(words[0]); // Name of boundary
  superPatch.setNbPatchesInSuperPatch(1); // Number of patches in each super patch is 1
  superPatch.getPatchIDs()[0] = patchCode;

  _superPatch.push_back(superPatch);

  patch.setPatchCode(patchCode);

  switch(bcType)
  {
   // Boundary condition type - Node or Element/Cell
   case 0: // Node type
           patch.setNbFacesInPatch(nbDataRecords-1);

     if (_dimension == 2)
     {
      // It is considered here that each couple of two following nodes in
      // Gambit file boundary condition build one face on boundary - 2D case
            getGambitWordsFromLine(fin,line,lineNb,words);
      firstFaceNode = (StringOps::from_str<CFuint>(words[0]))-1;
      for (CFuint i = 0; i < patch.getNbFacesInPatch(); ++i) {
              patch.getFaceData()[i].setCellID(i);
              patch.getFaceData()[i].setNbNodesInFace(2);

              getGambitWordsFromLine(fin,line,lineNb,words);
        secondFaceNode = (StringOps::from_str<CFuint>(words[0]))-1;

              patch.getFaceData()[i].getFaceNodes()[0] = firstFaceNode;
              patch.getFaceData()[i].getFaceNodes()[1] = secondFaceNode;

        firstFaceNode = secondFaceNode;
            }
           }
     else callGambitFileError("This type of boundary condition is treared just for 2D case.",lineNb,filepath);
           break;
   case 1: // Element/Cell type
     patch.setNbFacesInPatch(nbDataRecords);

     for (CFuint i = 0; i < patch.getNbFacesInPatch(); ++i) {
             getGambitWordsFromLine(fin,line,lineNb,words);

       cellID = StringOps::from_str<CFuint>(words[0]);
       gambitElemType = StringOps::from_str<CFuint>(words[1]);
       gambitFaceID = StringOps::from_str<CFuint>(words[2]);

       elemTypeOrder = getOrderChangeMap()(cellID-1,0); // Order of element type
       lineInElemTypeTable = getOrderChangeMap()(cellID-1,1); // Line in table of elements of element type
             nbNodesPerElemType = _elementType[elemTypeOrder].getNbNodesPerCell(); // Number of nodes of element type

       gambitElemTypeID = make_pair(gambitElemType,nbNodesPerElemType); // Gambit element type ID pair (Gambit element type, number of nodes per cell)
             faceDefinionTable = _mapFaces.find(gambitElemTypeID)->second; // Pointer to table of faces of element type

       nbFaceNodes = (*faceDefinionTable).nbCols(gambitFaceID-1); // Number of face nodes

       //calculating global order of elemnt in list of elements
             globalElemID = 0;
             for(CFuint j=0; j < elemTypeOrder; ++j)
       {
         globalElemID += _elementType[j].getNbCellsPerType();
       }
             globalElemID += lineInElemTypeTable; //global order of elemnt in list of elements

       patch.getFaceData()[i].setCellID(globalElemID);
             patch.getFaceData()[i].setNbNodesInFace(nbFaceNodes);

       for (CFuint j = 0; j < nbFaceNodes; ++j) {
         nodeOrderInElement = (*faceDefinionTable)(gambitFaceID-1,j); // Number of node corresponding to Gambit numbering of element nodes
               gambitFileNodeNumber = _elementType[elemTypeOrder].getNodeID(lineInElemTypeTable,nodeOrderInElement); // Number of node in Gambit file
         patch.getFaceData()[i].getFaceNodes()[j] = gambitFileNodeNumber;
             }
           }
           break;
   default: callGambitFileError("Wrong Gambit boundary condition format.",lineNb,filepath);
  }

  _patch.push_back(patch);

  // Checks if the ENDSECTION is reached
  getGambitWordsFromLine(fin,line,lineNb,words);
  if(words[0] != "ENDOFSECTION") callGambitFileError("Wrong number of boundary Elements/Nodes. ENDOFSECTION not reached.",lineNb,filepath);

  CFLogDebugMin(  "Gambit2CFmeshConverter::readBC() => end\n" );
}

////////////////////////////////////////////////////////////////////////////// 8

void Gambit2CFmeshConverter::readGambitFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
  path meshFile = change_extension(filepath, getOriginExtension());

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(meshFile);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;
  bool repeat_while = false;

  /**
  * Reads parameters and stores them in atributes of class Gambit2CFmeshConverter
  */
  readParameters(fin, meshFile, lineNb);

  // Goes to the header ENDOFSECTION - PUT IT INTO METHOD
  do
  {
    getGambitWordsFromLine(fin,line,lineNb,words);
    if(words.size()>0)
    {
      if(((words[0] != "ENDOFSECTION") && !fin.eof()) || (words[0] == "#")) repeat_while = true;
      else repeat_while = false;
    }
    else repeat_while = true;
  }
  while(repeat_while);

  // Read next line
  getGambitWordsFromLine(fin,line,lineNb,words);

  /**
  * Reads node coordinates
  */
  if(words[0] == "NODAL") readNodes(fin, meshFile, lineNb);
  else callGambitFileError("Wrong order of sections. NODAL COORDINATES section not found.",lineNb,meshFile);

  // Read next line
  getGambitWordsFromLine(fin,line,lineNb,words);

  /**
  * Reads element connectivity and stores it in atributes of class Gambit2CFmeshConverter
  */
  if(words[0] == "ELEMENTS/CELLS") readConnectivity(fin, meshFile, lineNb);
  else callGambitFileError("Wrong order of sections. ELEMENTS/CELLS section not found.",lineNb,meshFile);

  /**
  * Reads Groups and stores them
  */
  readGroups(fin, meshFile, lineNb);

  /// Only if boundary conditions are present
  for (CFuint i = 0; i < _nbSuperPatches; ++i )
  {
    // Read next line
    getGambitWordsFromLine(fin,line,lineNb,words);

    /**
    * Reads boundary conditions and stores them in atributes of class Gambit2CFmeshConverter
    */

    CFLogDebugMin(  "SuperPatch number = " << i << " => start\n" );
    if(words[0] == "BOUNDARY") readBC(fin, meshFile, lineNb, i);
    else callGambitFileError("Wrong order of sections. BOUNDARY CONDITIONS section not found.",lineNb,meshFile);

    CFLogDebugMin(  "SuperPatch number = " << i << " => end\n" );
  }

  fhandle->close();
}

////////////////////////////////////////////////////////////////////////////// 9

void Gambit2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  if(!_isFileRead) {
    try {
      readGambitFile(filepath);
      _isFileRead = true;
    }
    catch (Common::Exception& e) {
      CFout << e.what() << "\n";
      throw;
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 10

void Gambit2CFmeshConverter::writeContinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  /// @todo fix this for parallel mesh reading
  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  fout << "!NB_NODES "
       << _nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  fout << "!NB_STATES "
       << _nbUpdatableStates
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << _nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  /// @todo only first order for now
  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER1 << "\n";

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(_elementType[k].getNbNodesPerCell(),
                                      CFPolyOrder::ORDER1,
                                      _dimension) << " ";

  }
  fout << "\n";

  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << _elementType[k].getNbCellsPerType() << " ";
  }
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << _elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout <<_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    CFuint nbCellsPerType = _elementType[k].getNbCellsPerType();
    CFuint nbNodesPerCell = _elementType[k].getNbNodesPerCell();
    CFuint nbStatesPerCell  = nbNodesPerCell;

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j) {
        fout << _elementType[k].getNodeID(i,j) <<" " ;
      }
      for (CFuint j = 0; j < nbStatesPerCell; ++j) {
        fout << _elementType[k].getNodeID(i,j) <<" " ;
      }
      fout << "\n";
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 11

void Gambit2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << _isWithSolution << "\n";

  if (_isWithSolution) {
    fout.precision(14);
    fout << *_variables;
  }
}

////////////////////////////////////////////////////////////////////////////// 11

void Gambit2CFmeshConverter::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_NODE " << "\n";
  for (CFuint k = 0; k < _nbUpdatableNodes; ++k) {
    for (CFuint j = 0; j < _dimension; ++j)
    {
      fout.precision(14);
      fout.setf(ios::scientific,ios::floatfield);
      fout << _coordinate[k*_dimension+j] << " ";
    }
    fout << "\n";
  }
}

////////////////////////////////////////////////////////////////////////////// 12

void Gambit2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  /// @todo fix this for parallel mesh reading
  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  fout << "!NB_NODES "
       << _nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  /// @todo how to define nbNotUpdatableStates for FVM??
  fout << "!NB_STATES "
       << _nbCells
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << _nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  /// @todo only first order for now

  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER0 << "\n";
  // this CFPolyOrder::ORDER0 can only be set if here we know that we are dealing
  // with CellCenterFEM

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(_elementType[k].getNbNodesPerCell(),
              CFPolyOrder::ORDER1,
              _dimension) << " ";
  }
  fout << "\n";

  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << _elementType[k].getNbCellsPerType() << " ";
  }
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << _elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << 1 << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";

  CFuint countElem = 0;
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    const CFuint nbCellsPerType = _elementType[k].getNbCellsPerType();
    const CFuint nbNodesPerCell = _elementType[k].getNbNodesPerCell();

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j) {
        fout << _elementType[k].getNodeID(i,j) << " " ;
      }
      fout << countElem << "\n"; // cellID == stateID
      ++countElem;
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 13

void Gambit2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << _isWithSolution << "\n";

  CFuint nbVariables = PhysicalModelStack::getActive()->getNbEq();

  std::valarray<CFreal> averageState(nbVariables);

  /// @todo the following needs to be tested
  if (_isWithSolution) {
    for (CFuint k = 0; k < getNbElementTypes(); ++k) {
      CFuint nbCellsPerType = _elementType[k].getNbCellsPerType();
      CFuint nbNodesPerCell = _elementType[k].getNbNodesPerCell();

      cf_assert(nbCellsPerType > 0);
      for (CFuint i = 0; i < nbCellsPerType; ++i) {
  averageState = 0.0;
  for (CFuint j = 0; j < nbNodesPerCell; ++j) {
    const CFuint nodeID = _elementType[k].getNodeID(i,j);
    const vector<CFreal>& nodalState = _variables->getRow(nodeID);
    for(CFuint v = 0; v < nbVariables; ++v) {
      averageState[v] += nodalState[v];
    }
  }
  averageState /= nbNodesPerCell; // compute the average state
  for (CFuint iVar = 0; iVar < nbVariables; ++iVar) {
    fout.precision(14);
    fout << averageState[iVar] << " ";
  }
  fout << "\n";
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 13

void Gambit2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  if(_nbGroups > 1){
    fout << "!NB_GROUPS " << _nbGroups << "\n";

    //inner domain is split into groups
    for (CFuint iGroup = 0; iGroup < _nbGroups; ++iGroup) {

      const CFuint nbElementsInGroup = _groupCellIDs[iGroup].size();
      const std::string groupName = _groupNames[iGroup];

      fout << "!GROUP_NAME " << groupName << "\n";
      fout << "!GROUP_ELEM_NB " << nbElementsInGroup << "\n";
      fout << "!GROUP_ELEM_LIST" << "\n";

      for (CFuint iElem = 0; iElem < nbElementsInGroup; ++iElem) {
        fout << _groupCellIDs[iGroup][iElem] << "\n";
      }
    }
  }


  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;

  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < _nbPatches; ++iPatch) {
    CFuint codeID = _patch[iPatch].getPatchCode();
    mapCP[codeID] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << _nbSuperPatches << " " << nbUnsetTRS << "\n";

  const CFuint nbIgnoreTRSs = m_ignoreTRSNames.size();
  const CFuint nbTRSs = getNbSuperPatches() - nbIgnoreTRSs;
  fout << "!NB_TRSs " << nbTRSs << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = _superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = _superPatch[iTRS].getSuperPatchName();

    bool writeTRS = true;
    for (CFuint iName = 0; iName < nbIgnoreTRSs; ++iName) {
      if (nameTRS == m_ignoreTRSNames[iName]) {
  writeTRS = false;
      }
    }

    if (writeTRS) {
      fout << "!TRS_NAME " << nameTRS << "\n";
      fout << "!NB_TRs "<< nbTRsInTRS << "\n";
      fout << "!NB_GEOM_ENTS ";

      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
  CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
  CFuint curPatch = mapCP.find(patchID)->second;
  fout << _patch[curPatch].getNbFacesInPatch() << " ";
      }
      fout << "\n";
      fout << "!GEOM_TYPE Face" << "\n";
      fout << "!LIST_GEOM_ENT" << "\n";

      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

  CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
  CFuint curPatch = mapCP.find(patchID)->second;
  CFuint nbFacesInPatch = _patch[curPatch].getNbFacesInPatch();

  for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
    {

      CFuint nbNodesPerFace = _patch[curPatch].getFaceData()[iFace].getNbNodesInFace();
      CFuint nbStatesPerFace = nbNodesPerFace;

      fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
        fout << _patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
      }
      for (CFuint iState = 0; iState < nbStatesPerFace; ++iState) {
        fout << _patch[curPatch].getFaceData()[iFace].getFaceNodes()[iState] << " " ;
      }
      fout << "\n";
    }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 14

void Gambit2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;
  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < _nbPatches; ++iPatch) {
    CFuint codeP = _patch[iPatch].getPatchCode();
    mapCP[codeP] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << _nbSuperPatches << " " << nbUnsetTRS << "\n";
  const CFuint nbIgnoreTRSs = m_ignoreTRSNames.size();
  const CFuint nbTRSs = getNbSuperPatches() - nbIgnoreTRSs;
  fout << "!NB_TRSs " << nbTRSs << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = _superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = _superPatch[iTRS].getSuperPatchName();

    bool writeTRS = true;
    for (CFuint iName = 0; iName < nbIgnoreTRSs; ++iName) {
      if (nameTRS == m_ignoreTRSNames[iName]) {
        writeTRS = false;
      }
    }

    if (writeTRS) {
      fout << "!TRS_NAME " << nameTRS << "\n";
      fout << "!NB_TRs "<< nbTRsInTRS << "\n";
      fout << "!NB_GEOM_ENTS ";

      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
  CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
  CFuint curPatch = mapCP.find(patchID)->second;
  fout << _patch[curPatch].getNbFacesInPatch() << " ";
      }
      fout << "\n";
      fout << "!GEOM_TYPE Face" << "\n";
      fout << "!LIST_GEOM_ENT" << "\n";

      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

  CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
  CFuint curPatch = mapCP.find(patchID)->second;
  CFuint nbFacesInPatch = _patch[curPatch].getNbFacesInPatch();

  for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace) {
    CFuint nbNodesPerFace = _patch[curPatch].getFaceData()[iFace].getNbNodesInFace();
    CFuint nbStatesPerFace = 1;
    fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
    for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
      fout << _patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
    }
    /// @todo StateID hardlinked with the CellID
      fout << _patch[curPatch].getFaceData()[iFace].getCellID() << "\n";
  }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 15

void Gambit2CFmeshConverter::convertBack(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  throw Common::NotImplementedException (FromHere(),"Converting back from Gambit not implemented yet... have fun programing it.\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace Gambit2CFmesh

  } // end of namespace IO

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
