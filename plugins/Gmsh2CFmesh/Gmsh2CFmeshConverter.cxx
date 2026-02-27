// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/Stopwatch.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/CFMap.hh"
#include "Gmsh2CFmesh/Gmsh2CFmeshConverter.hh"
#include "Gmsh2CFmesh/Gmsh2CFmesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Gmsh2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Gmsh2CFmeshConverter,
               MeshFormatConverter,
               Gmsh2CFmeshModule,
               1>
gmsh2CFmeshConverterProvider("Gmsh2CFmesh");

//////////////////////////////////////////////////////////////////////////////

Gmsh2CFmeshConverter::Gmsh2CFmeshConverter (const std::string& name)
: MeshFormatConverter(name),
  _fileFormatVersion(0),
  _dimension(0),
  _order(0),
  _nbCells(0),
  _nbUpdatableNodes(0),
  _nbFaces(0),
  _nbPatches(0),
  _elementType(0),
  _patch(0),
  _superPatch(0),
  _update(static_cast<CFuint>(0),static_cast<CFuint>(0)),
  _coordinate(CFNULL),
  _variables(CFNULL),
  _isWithSolution(false),
  _isFileRead(false),
  _nbUpdatableStates(_nbUpdatableNodes),
  _nodesPerElemTypeTable(31),
  _orderPerElemTypeTable(31),
  _dimPerElemTypeTable(31),
  _mapNodeIdxPerElemTypeTable(31)
{
  // Build the nbNodes per ElemTypeTable
  _nodesPerElemTypeTable[0]  = 2;  // line
  _nodesPerElemTypeTable[1]  = 3;  // triangle
  _nodesPerElemTypeTable[2]  = 4;  // quadrangle
  _nodesPerElemTypeTable[3]  = 4;  // tetrahedron
  _nodesPerElemTypeTable[4]  = 8;  // hexahedron
  _nodesPerElemTypeTable[5]  = 6;  // prism
  _nodesPerElemTypeTable[6]  = 5;  // pyramid
  _nodesPerElemTypeTable[7]  = 3;  // P2 line
  _nodesPerElemTypeTable[8]  = 6;  // P2 triangle
  _nodesPerElemTypeTable[9]  = 9;  // P2 quadrangle
  _nodesPerElemTypeTable[10] = 10; // P2 tetrahedron
  _nodesPerElemTypeTable[11] = 27; // P2 hexahedron
  _nodesPerElemTypeTable[12] = 18; // P2 prism
  _nodesPerElemTypeTable[13] = 14; // P2 pyramid
  _nodesPerElemTypeTable[14] = 1;  // node
  _nodesPerElemTypeTable[15] = 8;  // P2 incomplete quadrangle
  _nodesPerElemTypeTable[16] = 20; // P2 incomplete hexahedron
  _nodesPerElemTypeTable[17] = 15; // P2 incomplete prism
  _nodesPerElemTypeTable[18] = 13; // P2 incomplete pyramid
  _nodesPerElemTypeTable[19] = 9;  // P3 incomplete triangle
  _nodesPerElemTypeTable[20] = 10; // P3 triangle
  _nodesPerElemTypeTable[21] = 12; // P4 incomplete triangle
  _nodesPerElemTypeTable[22] = 15; // P4 triangle
  _nodesPerElemTypeTable[23] = 15; // P5 incomplete triangle
  _nodesPerElemTypeTable[24] = 21; // P5 triangle
  _nodesPerElemTypeTable[25] = 4;  // P3 line
  _nodesPerElemTypeTable[26] = 5;  // P4 line
  _nodesPerElemTypeTable[27] = 6;  // P5 line
  _nodesPerElemTypeTable[28] = 20; // P3 tetrahedron
  _nodesPerElemTypeTable[29] = 35; // P4 tetrahedron
  _nodesPerElemTypeTable[30] = 56; // P5 tetrahedron

  // Build the order per ElemTypeTable
  _orderPerElemTypeTable[0]  = 1; // line
  _orderPerElemTypeTable[1]  = 1; // triangle
  _orderPerElemTypeTable[2]  = 1; // quadrangle
  _orderPerElemTypeTable[3]  = 1; // tetrahedron
  _orderPerElemTypeTable[4]  = 1; // hexahedron
  _orderPerElemTypeTable[5]  = 1; // prism
  _orderPerElemTypeTable[6]  = 1; // pyramid
  _orderPerElemTypeTable[7]  = 2; // P2 line
  _orderPerElemTypeTable[8]  = 2; // P2 triangle
  _orderPerElemTypeTable[9]  = 2; // P2 quadrangle
  _orderPerElemTypeTable[10] = 2; // P2 tetrahedron
  _orderPerElemTypeTable[11] = 2; // P2 hexahedron
  _orderPerElemTypeTable[12] = 2; // P2 prism
  _orderPerElemTypeTable[13] = 2; // P2 pyramid
  _orderPerElemTypeTable[14] = 0; // node
  _orderPerElemTypeTable[15] = 2; // P2 incomplete quadrangle
  _orderPerElemTypeTable[16] = 2; // P2 incomplete hexahedron
  _orderPerElemTypeTable[17] = 2; // P2 incomplete prism
  _orderPerElemTypeTable[18] = 2; // P2 incomplete pyramid
  _orderPerElemTypeTable[19] = 3; // P3 incomplete triangle
  _orderPerElemTypeTable[20] = 3; // P3 triangle
  _orderPerElemTypeTable[21] = 4; // P4 incomplete triangle
  _orderPerElemTypeTable[22] = 4; // P4 triangle
  _orderPerElemTypeTable[23] = 5; // P5 incomplete triangle
  _orderPerElemTypeTable[24] = 5; // P5 triangle
  _orderPerElemTypeTable[25] = 3; // P3 line
  _orderPerElemTypeTable[26] = 4; // P4 line
  _orderPerElemTypeTable[27] = 5; // P5 line
  _orderPerElemTypeTable[28] = 3; // P3 tetrahedron
  _orderPerElemTypeTable[29] = 4; // P4 tetrahedron
  _orderPerElemTypeTable[30] = 5; // P5 tetrahedron

  // Build the table which lists the dimension of each element
  _dimPerElemTypeTable[0]   = DIM_1D; // line
  _dimPerElemTypeTable[1]   = DIM_2D; // triangle
  _dimPerElemTypeTable[2]   = DIM_2D; // quadrangle
  _dimPerElemTypeTable[3]   = DIM_3D; // tetrahedron
  _dimPerElemTypeTable[4]   = DIM_3D; // hexahedron
  _dimPerElemTypeTable[5]   = DIM_3D; // prism
  _dimPerElemTypeTable[6]   = DIM_3D; // pyramid
  _dimPerElemTypeTable[7]   = DIM_1D; // P2 line
  _dimPerElemTypeTable[8]   = DIM_2D; // P2 triangle
  _dimPerElemTypeTable[9]   = DIM_2D; // P2 quadrangle
  _dimPerElemTypeTable[10]  = DIM_3D; // P2 tetrahedron
  _dimPerElemTypeTable[11]  = DIM_3D; // P2 hexahedron
  _dimPerElemTypeTable[12]  = DIM_3D; // P2 prism
  _dimPerElemTypeTable[13]  = DIM_3D; // P2 pyramid
  _dimPerElemTypeTable[14]  = DIM_0D; // node
  _dimPerElemTypeTable[15]  = DIM_2D; // P2 incomplete quadrangle
  _dimPerElemTypeTable[16]  = DIM_3D; // P2 incomplete hexahedron
  _dimPerElemTypeTable[17]  = DIM_3D; // P2 incomplete prism
  _dimPerElemTypeTable[18]  = DIM_3D; // P2 incomplete pyramid
  _dimPerElemTypeTable[19]  = DIM_2D; // P3 incomplete triangle
  _dimPerElemTypeTable[20]  = DIM_2D; // P3 triangle
  _dimPerElemTypeTable[21]  = DIM_2D; // P4 incomplete triangle
  _dimPerElemTypeTable[22]  = DIM_2D; // P4 triangle
  _dimPerElemTypeTable[23]  = DIM_2D; // P5 incomplete triangle
  _dimPerElemTypeTable[24]  = DIM_2D; // P5 triangle
  _dimPerElemTypeTable[25]  = DIM_1D; // P3 line
  _dimPerElemTypeTable[26]  = DIM_1D; // P4 line
  _dimPerElemTypeTable[27]  = DIM_1D; // P5 line
  _dimPerElemTypeTable[28]  = DIM_3D; // P3 tetrahedron
  _dimPerElemTypeTable[29]  = DIM_3D; // P4 tetrahedron
  _dimPerElemTypeTable[30]  = DIM_3D; // P5 tetrahedron

  // Build node indexes per ElemTypeTable
  for (CFuint iElemType = 0; iElemType < _nodesPerElemTypeTable.size(); ++iElemType)
  {
    _mapNodeIdxPerElemTypeTable[iElemType]
        .resize(_nodesPerElemTypeTable[iElemType]);
    for (CFuint iNode = 0; iNode < _nodesPerElemTypeTable[iElemType]; ++iNode)
    {
      _mapNodeIdxPerElemTypeTable[iElemType][iNode] = iNode;
    }
  }

  // correct the mapping of the node indices for some elements
  /// @note this should be checked for quadratic elements not listed below, as the node order might differ,
  /// and some elements are not yet implemented in COOLFluiD!
  /// okay for quadratic line, quadratic quadrangle
  /// @note it appears that the order the nodes in the face centers are written to the mesh file is not always the same...
  /// @note this is due to a bug in Gmsh (up to version 2.0.8). It will be fixed in a forthcoming release of Gmsh
  // quadratic hexahedron
  _mapNodeIdxPerElemTypeTable[11][ 8] =  8;
  _mapNodeIdxPerElemTypeTable[11][ 9] = 11;
  _mapNodeIdxPerElemTypeTable[11][10] = 13;
  _mapNodeIdxPerElemTypeTable[11][11] =  9;
  _mapNodeIdxPerElemTypeTable[11][12] = 20;
  _mapNodeIdxPerElemTypeTable[11][13] = 10;
  _mapNodeIdxPerElemTypeTable[11][14] = 21;
  _mapNodeIdxPerElemTypeTable[11][15] = 12;
  _mapNodeIdxPerElemTypeTable[11][16] = 23;
  _mapNodeIdxPerElemTypeTable[11][17] = 14;
  _mapNodeIdxPerElemTypeTable[11][18] = 24;
  _mapNodeIdxPerElemTypeTable[11][19] = 15;
  _mapNodeIdxPerElemTypeTable[11][20] = 22;
  _mapNodeIdxPerElemTypeTable[11][21] = 26;
  _mapNodeIdxPerElemTypeTable[11][22] = 16;
  _mapNodeIdxPerElemTypeTable[11][23] = 18;
  _mapNodeIdxPerElemTypeTable[11][24] = 19;
  _mapNodeIdxPerElemTypeTable[11][25] = 17;
  _mapNodeIdxPerElemTypeTable[11][26] = 25;
  
  // quadratic Prism
  _mapNodeIdxPerElemTypeTable[12][ 6] =  6;
  _mapNodeIdxPerElemTypeTable[12][ 7] =  9;
  _mapNodeIdxPerElemTypeTable[12][ 8] =  7;
  _mapNodeIdxPerElemTypeTable[12][ 9] =  8;
  _mapNodeIdxPerElemTypeTable[12][10] = 15;
  _mapNodeIdxPerElemTypeTable[12][11] = 10;
  _mapNodeIdxPerElemTypeTable[12][12] = 17;
  _mapNodeIdxPerElemTypeTable[12][13] = 11;
  _mapNodeIdxPerElemTypeTable[12][14] = 16;
  _mapNodeIdxPerElemTypeTable[12][15] = 12;
  _mapNodeIdxPerElemTypeTable[12][16] = 14;
  _mapNodeIdxPerElemTypeTable[12][17] = 13;



  _mapNodeIdxPerElemTypeTable[10][ 7] =  9;
  _mapNodeIdxPerElemTypeTable[10][ 9] = 7;

  // incomplete quadratic hexahedron
  _mapNodeIdxPerElemTypeTable[16][ 8] =  8;
  _mapNodeIdxPerElemTypeTable[16][ 9] = 11;
  _mapNodeIdxPerElemTypeTable[16][10] = 13;
  _mapNodeIdxPerElemTypeTable[16][11] =  9;
  _mapNodeIdxPerElemTypeTable[16][12] = 10;
  _mapNodeIdxPerElemTypeTable[16][13] = 12;
  _mapNodeIdxPerElemTypeTable[16][14] = 14;
  _mapNodeIdxPerElemTypeTable[16][15] = 15;
  _mapNodeIdxPerElemTypeTable[16][16] = 16;
  _mapNodeIdxPerElemTypeTable[16][17] = 18;
  _mapNodeIdxPerElemTypeTable[16][18] = 19;
  _mapNodeIdxPerElemTypeTable[16][19] = 17;
}

//////////////////////////////////////////////////////////////////////////////

Gmsh2CFmeshConverter::~Gmsh2CFmeshConverter()
{
  deletePtr(_coordinate);
  deletePtr(_variables);
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::checkFormat(const boost::filesystem::path& filepath)
{

  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path meshFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
  path meshFile = change_extension(filepath, getOriginExtension());
#endif

  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& fin = (*fhandle)->open(meshFile);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;

//  CFuint elemNb;
  CFuint elemType;
//  CFuint phyRegion;
//  CFuint elemRegion;

  _dimension = DIM_0D;

  // check
  getGmshWordsFromLine(fin,line,lineNb,words);

  // For file formats, see www.geuz.org/gmsh
  // Old file format:
  // $NOD
  // number-of-nodes
  // $ENDNOD
  // $ELM
  // number-of-elements
  // elm-number elm-type reg-phys re-elem number-of-nodes node-number-list
  // ...
  // $ENDELM

  if ( (words.size() == 1) && (words[0] == "$NOD") )
  {
    CFout << "The file seems to have Gmsh version 1 format ...\n";

    getGmshWordsFromLine(fin,line,lineNb,words);
    _nbUpdatableNodes = StringOps::from_str<CFint>(words[0]);
    for(CFuint i = 0; i < _nbUpdatableNodes; ++i)
    {
      getline(fin,line);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if(words[0] != "$ENDNOD")
    {
      callGmshFileError("Malformed file format: $ENDNOD statement missing",
                                                 lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);

    if(words[0] != "$ELM")
    {
      callGmshFileError("Malformed file format: $ELM statement missing",
                                                 lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    const CFuint nbElements = StringOps::from_str<CFint>(words[0]);

    for(CFuint i = 0; i < nbElements; ++i)
    {
      getGmshWordsFromLine(fin,line,lineNb,words);
      elemType = StringOps::from_str<CFuint>(words[1]);
      _dimension = std::max(_dimension,_dimPerElemTypeTable[elemType-1]);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if(words[0] != "$ENDELM")
    {
      callGmshFileError("Malformed file format: $ENDELM statement missing",
                                                 lineNb,meshFile);
    }

    _fileFormatVersion = 1;

  } // Check of file format version 1

  // Remark: acoording to the manual, the order of the sections in the gmsh file
  // is not fixed, but for the moment, we suppose that the sections
  // are organized as follows:
  // $MeshFormat
  // ... info about mesh format
  // $EndMeshFormat
  // $PhysicalNames
  // ... list of physical entities
  // $EndPhysicalNames
  // $Nodes
  // ... list of mesh nodes
  // $EndNodes
  // $Elements
  // ... list of all elements in the mesh
  // $EndElements

  else if ( (words.size() == 1) && (words[0] == "$MeshFormat") )
  {
    getGmshWordsFromLine(fin,line,lineNb,words);
    const CFreal versionNumber = StringOps::from_str<CFreal>(words[0]);

    if (versionNumber >= 2.0)
    {
      CFout << "The file seems to have Gmsh version " << versionNumber << " format\n";
    }
    getGmshWordsFromLine(fin,line,lineNb,words);
    if (words[0] != "$EndMeshFormat")
    {
      callGmshFileError("Malformed file format: $EndMeshFormat statement missing",
                                                 lineNb,meshFile);
    }
    getGmshWordsFromLine(fin,line,lineNb,words);

    if (words[0] != "$PhysicalNames")
    {
      callGmshFileError("Malformed file format: $PhysicalNames statement missing",
                                                lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);

    const CFuint nbPhysicalNames = StringOps::from_str<CFuint>(words[0]);

    for(CFuint i = 0; i < nbPhysicalNames; ++i)
    {
      getline(fin,line);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if (words[0] != "$EndPhysicalNames")
    {
      callGmshFileError("Malformed file format: $EndPhysicalNames statement missing",
                                                lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if (words[0] != "$Nodes")
    {
      callGmshFileError("Malformed file format: $Nodes statement missing",
                                                lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    const CFuint nbNodes = StringOps::from_str<CFuint>(words[0]);

    for(CFuint i = 0; i < nbNodes; ++i)
    {
      getline(fin,line);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if (words[0] != "$EndNodes")
    {
      callGmshFileError("Malformed file format: $EndNodes statement missing",
                                                lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if (words[0] != "$Elements")
    {
      callGmshFileError("Malformed file format: $Elements statement missing",
                                                lineNb,meshFile);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    const CFuint nbElements = StringOps::from_str<CFuint>(words[0]);

    for(CFuint i = 0; i < nbElements; ++i)
    {
      getGmshWordsFromLine(fin,line,lineNb,words);
      elemType = StringOps::from_str<CFuint>(words[1]);
      _dimension = std::max(_dimension,_dimPerElemTypeTable[elemType-1]);
    }

    getGmshWordsFromLine(fin,line,lineNb,words);
    if (words[0] != "$EndElements")
    {
      callGmshFileError("Malformed file format: $EndElements statement missing",
                                                lineNb,meshFile);
    }

    _fileFormatVersion = 2;

  } // Check of file format version 2

  // check nb Elements,
  // nb nodes, nb Boundary Faces
  // nb patches

  ///@todo add many more checks

  (*fhandle)->close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::callGmshFileError(const std::string& msg,
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

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::getGmshWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::readGmshFileVersion1(const boost::filesystem::path& filepath)
{
   CFAUTOTRACE;

   using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
   path meshFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
   path meshFile = change_extension(filepath, getOriginExtension());
#endif

   Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
   ifstream& fin = (*fhandle)->open(meshFile);

   CFuint lineNb = 0;
   std::string line;
   vector<std::string> words;

   // read general information
   // check
   getGmshWordsFromLine(fin,line,lineNb,words);
   if(words.size() != 1)
      callGmshFileError("Wrong number of parameters.",lineNb,meshFile);

   if(words[0] != "$NOD")
      callGmshFileError("Malformed file format",lineNb,meshFile);

   // Read Number of nodes
   getGmshWordsFromLine(fin,line,lineNb,words);
   _nbUpdatableNodes = StringOps::from_str<CFint>(words[0]);

   // Get dimension from first node
   ///@todo Gmsh always stores 3D... 2D is with Z=0
   //_dimension = 3;
   //Update: _dimension is already set by checkFormat(...)

   /// A map is needed between the Gmsh numbering and an ordered numbering!!
   /// The numbering of the nodes in Gmsh is ARBITRARY, i.e. the nodes
   /// in the mesh file are not necessarily listed in order 1,2,..,n

   CFMap<CFuint,CFuint> gmshNodesNumberingMap(_nbUpdatableNodes);
   // unused //  CFuint gmshNumbering;

   _coordinate = new Table<CFreal>(_nbUpdatableNodes, _dimension);

   for (CFuint j = 0; j < _nbUpdatableNodes; ++j)
   {
      getGmshWordsFromLine(fin,line,lineNb,words);
      // insert map between the node number and the ordered node number
      gmshNodesNumberingMap.insert(StringOps::from_str<CFuint>(words[0]),j);
      for (CFuint i = 0; i < _dimension; ++i)
      {
         (*_coordinate)(j,i) = StringOps::from_str<CFreal>(words[i+1]);
      }
   }

   gmshNodesNumberingMap.sortKeys();
   // check
   getGmshWordsFromLine(fin,line,lineNb,words);
   if(words.size() != 1) callGmshFileError("Wrong number of parameters.",lineNb,meshFile);
   if(words[0] != "$ENDNOD") callGmshFileError("Malformed file format, there should be $ENDNOD at line ",lineNb,meshFile);
   // check
   getGmshWordsFromLine(fin,line,lineNb,words);
   if(words.size() != 1) callGmshFileError("Wrong number of parameters.",lineNb,meshFile);
   if(words[0] != "$ELM") callGmshFileError("Malformed file format, waiting for $ELM descriptor at line ",lineNb,meshFile);
   // Read Number of Elements
   getGmshWordsFromLine(fin,line,lineNb,words);
   if(words.size() != 1) callGmshFileError("Wrong number of parameters.",lineNb,meshFile);
   _nbCells = StringOps::from_str<CFint>(words[0]);
   // Elements are not grouped by element types...
   // First loop over all the elements to know the number of element types
   // Also get the number of Physical Regions (<->SP) and elemnt Regions (<->Patches)
   const CFuint maxNbElementTypes = _nodesPerElemTypeTable.size();
   vector<CFuint> nbElemPerType(maxNbElementTypes,0);

   vector<CFuint> phyReg(0);
   vector<CFuint> elemReg(0);
   vector<CFuint> elemsPerPhyReg(0);
   vector<CFuint> elemsPerElemReg(0);
   // Build some temporaries
   CFuint nbNodes;
   CFuint elemNb;
   CFuint elemType;
   CFuint phyRegion;
   CFuint elemRegion;
   // unused // CFuint nbNodesInCell;
   bool found = false;
   // Temporary int to count the number of elements in InnerCells
   CFuint totalNbElements = 0;
   for (CFuint i = 0; i < _nbCells; ++i){
      getGmshWordsFromLine(fin,line,lineNb,words);
      elemNb      = StringOps::from_str<CFint>(words[0]);
      elemType    = StringOps::from_str<CFint>(words[1]);
      phyRegion   = StringOps::from_str<CFint>(words[2]);
      elemRegion  = StringOps::from_str<CFint>(words[3]);
 
      bool isPoint(elemType == 15);
      bool isEdge(false);
      if((_dimension == DIM_3D) && (elemType == 8)) isEdge = true;

      if(!isPoint && !isEdge)
      {
         if (phyRegion == _inFieldSP){
            nbElemPerType[elemType-1] += 1;
            totalNbElements += 1;
         }

         /// Count the number of PhysRegions and NbElements per Region
         // Check if phyRegion already exists.
         found = false;
         for (CFuint i=0; i < phyReg.size(); ++i){
            if (phyReg[i] == phyRegion){
               found = true;
               elemsPerPhyReg[i] += 1;
            }
         }
         // If new phyRegion, add it to the list.
         if (!found){
            phyReg.push_back(phyRegion);
            elemsPerPhyReg.push_back(1);
         }

         /// Count the number of ElementRegions and NbElements per Region
         // Check if ElemRegion already exists.
         found = false;
         for (CFuint i=0; i < elemReg.size(); ++i){
            if (elemReg[i] == elemRegion){
               found = true;
               elemsPerElemReg[i] += 1;
            }
         }
         // If new ElemRegion, add it to the list.
         if (!found){
            elemReg.push_back(elemRegion);
            elemsPerElemReg.push_back(1);
         }
      }
     
   }
   // Now we know the number of:
   // - Elements
   // - ElementTypes and Nb of Elements per Type
   // - Regions and Nb of Elements per Region
   /// This allows the allocation of the memory
   // Count the number of elementTypes with non-zero elements
   CFuint nbElementTypes = 0;
   for (CFuint i=0; i < maxNbElementTypes; ++i){
      if (nbElemPerType[i] > 0) ++nbElementTypes;
   }
   // Read and allocate element type information
   _elementType.resize(nbElementTypes);
   // For each elementType present in the mesh, set the typeID
   CFuint k=0;
   for (CFuint i = 0; i < maxNbElementTypes; ++i){
      if (nbElemPerType[i] > 0){
         _elementType[k].setTypeID(i);
         ++k;
      }
   }
//
   // For each elementType present in the mesh, set the nbNodesPerCell, the order and NbCellsPerType
   for (CFuint k = 0; k < getNbElementTypes(); ++k) {
      _elementType[k].setNbNodesPerCell(_nodesPerElemTypeTable[_elementType[k].getTypeID()]);
      _elementType[k].setOrderPerCell  (_orderPerElemTypeTable[_elementType[k].getTypeID()]);
      _order = _order > _orderPerElemTypeTable[_elementType[k].getTypeID()] ?
               _order : _orderPerElemTypeTable[_elementType[k].getTypeID()];
      _elementType[k].setNbCellsPerType(nbElemPerType[_elementType[k].getTypeID()]);
      _elementType[k].setCurrentIndex(0);
   }
   // Create storage space for the connectivity tables
   for (CFuint k = 0; k < nbElementTypes; ++k) {
      _elementType[k].createCellNodeConnectivity();
   }
   // allocate patch storage space
   _nbPatches = elemReg.size();
   _patch.resize(_nbPatches);
   for (CFuint k = 0; k < _nbPatches; ++k) {
      _patch[k].setPatchCode(elemReg[k]);
      _patch[k].setNbFacesInPatch(elemsPerElemReg[k]);
      _patch[k].setCurrentIndex(0);
   }
   /// Re-loop over the elements to effectively read the elements
   //Change the line Number to reloop over the elements
   lineNb -= _nbCells;
   /// Go back to the beginning of the elements
   /// @todo try to do this without closing/opening...
   (*fhandle)->close();
   (*fhandle)->open(meshFile);
   for (CFuint i = 0; i < lineNb; ++i){
      getline(fin,line);
   }
   CFuint elemTypeIdx = 0;
   CFuint inTypeIdx;
   CFuint inPatchIdx;
   CFuint patchIdx = 0;
   CFuint localNumbering;
   for (CFuint i = 0; i < _nbCells; ++i){

      getGmshWordsFromLine(fin,line,lineNb,words);

      elemNb      = StringOps::from_str<CFint>(words[0]);
      elemType    = StringOps::from_str<CFint>(words[1]);
      phyRegion   = StringOps::from_str<CFint>(words[2]);
      elemRegion  = StringOps::from_str<CFint>(words[3]);
      nbNodes     = StringOps::from_str<CFint>(words[4]);
      CFuint ordre = _orderPerElemTypeTable[elemType-1];  
      bool isPoint(elemType == 15);
      bool isEdge(false);
      if((_dimension == DIM_3D) && (ordre== 2) && (elemType ==8)) isEdge = true;
	if((_dimension == DIM_3D) && (ordre== 1) && (elemType ==1)) isEdge = true;	
      if(!isPoint && !isEdge)
      {
         ///Only store the InField into elemType
         if (phyRegion == _inFieldSP){
            //What is the type of the cell??
            found = false;
            for (CFuint j=0; j < _elementType.size(); ++j){
               if (_elementType[j].getTypeID() == elemType-1){
                  elemTypeIdx = j;
                  found = true;
               }
            }

            inTypeIdx = _elementType[elemTypeIdx].getCurrentIndex();

            for (CFuint j = 0; j < nbNodes; ++j){
               localNumbering = gmshNodesNumberingMap.find(StringOps::from_str<CFint>(words[5+j]));
               _elementType[elemTypeIdx].getTableConnectivity()(inTypeIdx,j) = localNumbering;
            }
            _elementType[elemTypeIdx].setCurrentIndex(inTypeIdx+1);
         }

         if (phyRegion != _inFieldSP){
            // read the Element Region information
            // What is the vector index of the patch
            for (CFuint j=0; j < _patch.size(); ++j){
               if (_patch[j].getPatchCode() == elemRegion){
                  patchIdx = j;
               }
            }
            inPatchIdx = _patch[patchIdx].getCurrentIndex();
            /*
             * cellID is never defined !!!
             * Instead we create a node->elementmap that will be used
             * when writing the TR entities to find out the cell number
             */

            // Set to zero to avoid undefined values
            _patch[patchIdx].getFaceData()[inPatchIdx].setCellID(0);

            _patch[patchIdx].getFaceData()[inPatchIdx].setNbNodesInFace(nbNodes);
            for (CFuint j = 0; j < nbNodes; ++j) {
               localNumbering = gmshNodesNumberingMap.find(StringOps::from_str<CFint>(words[5+j]));
               _patch[patchIdx].getFaceData()[inPatchIdx].getFaceNodes()[j] = localNumbering;
            }
            _patch[patchIdx].setCurrentIndex(inPatchIdx + 1);

         }
      }
   }
   _nbCells = totalNbElements;
  (*fhandle)->close();
  delete fhandle;
}


//////////////////////////////////////////////////////////////////////////////


void Gmsh2CFmeshConverter::readGmshFileVersion2(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path meshFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
  path meshFile = change_extension(filepath, getOriginExtension());
#endif

  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& fin = (*fhandle)->open(meshFile);

  CFuint lineNb = 0;
  std::string line = "";
  vector<std::string> words;
  map<CFuint,CFuint>      physTagToDim;
  map<CFuint,std::string> physTagToName;

  while(line != "$PhysicalNames") getline(fin,line);
  getGmshWordsFromLine(fin,line,lineNb,words);

  const CFuint nbPhysicalNames = StringOps::from_str<CFuint>(words[0]);


  // NOTE: we suppose that the physical names section has the following format:
  // nb_of_phys_regions
  // dim index_of_phys_region region_name
  // dim index_of_phys_region region_name
  // ...
  // and the index_of_phys_region is a number such that 1 <= index_of_phys_region <= nb_of_phys_regions, i.e. the 
  // numbering starts by 1!
  for(CFuint i = 0; i < nbPhysicalNames; ++i)
  {
    getGmshWordsFromLine(fin,line,lineNb,words);
    physTagToDim.insert(make_pair<CFuint,CFuint>( StringOps::from_str<CFuint>(words[1]),StringOps::from_str<CFuint>(words[0])) );

    // Clip off the leading and trailing quote
    std::string regionName = words[2].substr(1,words[2].length()-2);
    physTagToName[StringOps::from_str<CFuint>(words[1])] = regionName;
  }

  getGmshWordsFromLine(fin,line,lineNb,words);
  if (words[0] != "$EndPhysicalNames")
  {
    callGmshFileError("Malformed file format: $EndPhysicalNames statement missing",
                                                lineNb,meshFile);
  }

  while(line != "$Nodes") getline(fin,line);

  getGmshWordsFromLine(fin,line,lineNb,words);

  _nbUpdatableNodes = StringOps::from_str<CFint>(words[0]);

  CFMap<CFuint,CFuint> gmshNodesNumberingMap(_nbUpdatableNodes);
   
  _coordinate = new Table<CFreal>(_nbUpdatableNodes, _dimension);

  for (CFuint j = 0; j < _nbUpdatableNodes; ++j)
  {
    getGmshWordsFromLine(fin,line,lineNb,words);
    // insert map between the node number and the ordered node number
    gmshNodesNumberingMap.insert(StringOps::from_str<CFuint>(words[0]),j);
    for (CFuint i = 0; i < _dimension; ++i)
    {
      (*_coordinate)(j,i) = StringOps::from_str<CFreal>(words[i+1]);
    }
  }

  gmshNodesNumberingMap.sortKeys();
  
  line = "";

  while(line != "$Elements") getline(fin,line);

  getGmshWordsFromLine(fin,line,lineNb,words);

  const CFuint totalNbElements = StringOps::from_str<CFint>(words[0]);

   // Elements are not grouped by element types...
   // First loop over all the elements to know the number of element types
   // Also get the number of physical regions (<->SP)  and elementary regions (<->Patches)
   const CFuint maxNbElementTypes = _nodesPerElemTypeTable.size();
   vector<CFuint> nbElemPerType(maxNbElementTypes,0);
  
   const CFuint nbRegions = physTagToDim.size();
 
   //vector<CFuint> phyReg(nbRegions,0);
   //vector<CFuint> elemsPerPhyReg(nbRegions,0);
   //vector<CFuint> elemReg(0);
   //vector<CFuint> elemsPerElemReg(0);
   //
   
   //The first number in each pair is the physical region index,
   //the second is the number of elements in this physical region
   map<CFuint,CFuint> elemsPerPhyReg;
   //The first number in each pair is the elementary region index,
   //the second is the number of elements in this region
   //map<CFuint,CFuint> elemsPerElemReg;

   // Build some temporaries
   CFuint nbNodes;
   CFuint elemNb;
   CFuint elemType;
   CFuint nbTags;
   CFuint phyRegion;
   CFuint elemRegion;
  
  // Temporary int to count the number of elements in InnerCells
  _nbCells = 0;
  for (CFuint i = 0; i < totalNbElements; ++i)
  {
    getGmshWordsFromLine(fin,line,lineNb,words);
    elemNb      = StringOps::from_str<CFuint>(words[0]);
    elemType    = StringOps::from_str<CFuint>(words[1]);
    nbTags      = StringOps::from_str<CFuint>(words[2]);
    phyRegion   = StringOps::from_str<CFuint>(words[3]);

    if (nbTags > 1)
    {
      elemRegion  = StringOps::from_str<CFuint>(words[4]);
    }
    else
    {
      callGmshFileError("The element connectivity has missing elementary tags",lineNb,meshFile);
    }

    const map<CFuint,CFuint>::const_iterator regionIt = physTagToDim.find(phyRegion);

    if ( regionIt == physTagToDim.end() )
    {
      callGmshFileError("Element with unknown physical region number",lineNb,meshFile);
    }

    elemsPerPhyReg[regionIt->first] += 1;

    // Count volume cells
    if ( _dimPerElemTypeTable[elemType-1] == _dimension )
    {
      _nbCells++;
      nbElemPerType[elemType-1] += 1;
    }

  } // Loop over all elements in the mesh

  // Count the number of elementTypes with non-zero elements
  CFuint nbElementTypes = 0;
  for (CFuint i=0; i < maxNbElementTypes; ++i)
  {
    if (nbElemPerType[i] > 0) ++nbElementTypes;
  }

  // Now we know the number of:
  // - Elements
  // - ElementTypes and Nb of Elements per Type
  // - Regions and Nb of Elements per Region
  // This allows the allocation of the memory
   
   // Read and allocate element type information
   _elementType.resize(nbElementTypes);
   // For each elementType present in the mesh, set the typeID
   CFuint k=0;
   for (CFuint i = 0; i < maxNbElementTypes; ++i)
   {
     if (nbElemPerType[i] > 0)
     {
        _elementType[k].setTypeID(i);
        ++k;
     }
   }

   // For each elementType present in the mesh, set the nbNodesPerCell, the order and NbCellsPerType
   for (CFuint k = 0; k < getNbElementTypes(); ++k)
   {
     _elementType[k].setNbNodesPerCell(_nodesPerElemTypeTable[_elementType[k].getTypeID()]);
     _elementType[k].setOrderPerCell  (_orderPerElemTypeTable[_elementType[k].getTypeID()]);
     _order = _order > _orderPerElemTypeTable[_elementType[k].getTypeID()] ?  _order : _orderPerElemTypeTable[_elementType[k].getTypeID()];
     _elementType[k].setNbCellsPerType(nbElemPerType[_elementType[k].getTypeID()]);
     _elementType[k].setCurrentIndex(0);
   }

   // Create storage space for the connectivity tables
   for (CFuint k = 0; k < nbElementTypes; ++k) 
   {
     _elementType[k].createCellNodeConnectivity();
   }


   // allocate patch storage space
   _nbPatches = 0;

    for(map<CFuint,CFuint>::const_iterator regionIt = physTagToDim.begin(); regionIt != physTagToDim.end(); ++regionIt)
    {
      if ( regionIt->second == (_dimension-1) )
      {
        _nbPatches++;
      }
    }

   _patch.resize(_nbPatches);
   _superPatch.resize(_nbPatches);
   k = 0;
   for (map<CFuint,CFuint>::const_iterator regionIt = physTagToDim.begin(); regionIt != physTagToDim.end(); ++regionIt)
   {
      if ( regionIt->second == (_dimension-1) )
      {
        _patch[k].setPatchCode(regionIt->first);

        const map<CFuint,CFuint>::const_iterator patchIt = elemsPerPhyReg.find(regionIt->first);
        _patch[k].setNbFacesInPatch(patchIt->second);
        _patch[k].setCurrentIndex(0);

        const map<CFuint,string>::const_iterator superPatchNameIt = physTagToName.find(regionIt->first);

        _superPatch[k].setSuperPatchName(superPatchNameIt->second);
        // For the moment, we'll suppose that each "superpatch" has only one boundary patch
        // in gmsh file format ver. 2
        _superPatch[k].setNbPatchesInSuperPatch(1);

        std::valarray<CFuint>& patchIDs = _superPatch[k].getPatchIDs();
        patchIDs.resize(1);
        patchIDs[0] = regionIt->first;

        k++;
      }
   }

   /// Re-loop over the elements to effectively read the elements
   /// Go back to the beginning of the elements
   /// @todo try to do this without closing/opening...
   (*fhandle)->close();
   (*fhandle)->open(meshFile);

   line = "";
   while( line != "$Elements" )
   {
     getline(fin,line);
   }
   getline(fin,line);

   CFuint elemTypeIdx = 0;
   CFuint inTypeIdx;
   CFuint inPatchIdx;
   CFuint patchIdx = 0;
   CFuint localNumbering;

   _nbCells = 0;
   for (CFuint i = 0; i < totalNbElements; ++i)
   {

     getGmshWordsFromLine(fin,line,lineNb,words);

     elemNb      = StringOps::from_str<CFint>(words[0]);
     elemType    = StringOps::from_str<CFint>(words[1]);
     nbTags      = StringOps::from_str<CFint>(words[2]);
     phyRegion   = StringOps::from_str<CFint>(words[3]);
//      elemRegion  = StringOps::from_str<CFint>(words[3]);
     nbNodes     = _nodesPerElemTypeTable[elemType-1];

     // If this is a volume cell (has the same dimension as the mesh)
     if ( _dimPerElemTypeTable[elemType-1] == _dimension )
     {
       _nbCells++;

       for (CFuint j=0; j < _elementType.size(); ++j)
       {
          if (_elementType[j].getTypeID() == elemType-1)
          {
             elemTypeIdx = j;
          }
       }

       inTypeIdx = _elementType[elemTypeIdx].getCurrentIndex();

       for (CFuint j = 0; j < nbNodes; ++j)
       {
          localNumbering = gmshNodesNumberingMap.find(StringOps::from_str<CFint>(words[3+nbTags+j]));
          _elementType[elemTypeIdx].getTableConnectivity()(inTypeIdx,j) = localNumbering;
       }
       _elementType[elemTypeIdx].setCurrentIndex(inTypeIdx+1);

     }

     else if ( _dimPerElemTypeTable[elemType-1] == (_dimension-1) )
     {
       // read the Element Region information
       // What is the vector index of the patch
       for (CFuint j=0; j < _patch.size(); ++j){
          if (_patch[j].getPatchCode() == phyRegion)
          {
             patchIdx = j;
          }
       }
       inPatchIdx = _patch[patchIdx].getCurrentIndex();
       /*
        * cellID is never defined !!!
        * Instead we create a node->elementmap that will be used
        * when writing the TR entities to find out the cell number
        */

       // Set to zero to avoid undefined values
       _patch[patchIdx].getFaceData()[inPatchIdx].setCellID(0);

       _patch[patchIdx].getFaceData()[inPatchIdx].setNbNodesInFace(nbNodes);
       for (CFuint j = 0; j < nbNodes; ++j)
       {
          localNumbering = gmshNodesNumberingMap.find(StringOps::from_str<CFint>(words[5+j]));
          _patch[patchIdx].getFaceData()[inPatchIdx].getFaceNodes()[j] = localNumbering;
       }
       _patch[patchIdx].setCurrentIndex(inPatchIdx + 1);

     }

  } // Loop over cells


  (*fhandle)->close();
  delete fhandle;	
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::readSPInnerCells(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
   path fileSP = boost::filesystem::path(filepath).replace_extension(".SP");
#else
  path fileSP = change_extension(filepath,".SP");
#endif
  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& fin = (*fhandle)->open(fileSP);

  std::string line;
  vector<std::string> words;

  getline(fin,line);
  words = Common::StringOps::getWords(line);

  if(words.size() != 1) {
    throw BadFormatException (FromHere(),"Bad number of parameters in file " + fileSP.string());
  }

  const CFuint nbSuperPatches = StringOps::from_str<CFuint>(words[0]);

  bool found(false);
  for (CFuint i = 0; i < nbSuperPatches && !found; ++i) {

    getline(fin,line);
    words = Common::StringOps::getWords(line);

    std::string superPatchName = words[0];
    if (superPatchName == "InField") {
      found = true;
      getline(fin,line);
      words = Common::StringOps::getWords(line);
      _inFieldSP = StringOps::from_str<CFint>(words[0]);
      }
    else getline(fin,line);
    }

  if(!found) throw BadFormatException (FromHere(),"No InField defined in SP file");

  (*fhandle)->close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::readSPFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path fileSP = boost::filesystem::path(filepath).replace_extension(".SP");
#else
  path fileSP = change_extension(filepath,".SP");
#endif

  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& fin = (*fhandle)->open(fileSP);

  std::string line;
  vector<std::string> words;

  getline(fin,line);
  words = Common::StringOps::getWords(line);

  if(words.size() != 1) {
    throw BadFormatException (FromHere(),"Bad number of parameters in file " + fileSP.string());
  }

  CFuint nbSuperPatches = StringOps::from_str<CFuint>(words[0]);

  if(nbSuperPatches > _nbPatches) {
    throw BadFormatException (FromHere(),"Number of Physical Regions bigger than Element Regions " + fileSP.string());
  }

  /// Here we remove the InField from the list of SuperPatches
  _superPatch.resize(nbSuperPatches-1);

  // unused //  bool found(false);
  for (CFuint i = 0; i < nbSuperPatches-1; ++i) {

    getline(fin,line);
    words = Common::StringOps::getWords(line);

    std::string superPatchName = words[0];
    CFint nbPatchesInSuperPatch = 0;

    if(superPatchName != "InField") {
      if(words.size() != 2) {
        throw BadFormatException (FromHere(),"Bad number of parameters in file " + fileSP.string());
        }
      nbPatchesInSuperPatch = StringOps::from_str<CFint>(words[1]);
      if(static_cast<CFuint>(nbPatchesInSuperPatch) > _nbPatches ||
         nbPatchesInSuperPatch <= 0) {
        throw BadFormatException (FromHere(),"Bad number of Patches: "
             + words[1]
             + " in SuperPatch "
             + fileSP.string());
        }
      }

    if (superPatchName != "InField") {

    _superPatch[i].setSuperPatchName(superPatchName);
    _superPatch[i].setNbPatchesInSuperPatch(nbPatchesInSuperPatch);

    CFuint foundIDs = 0;
    while (foundIDs < static_cast<CFuint>(nbPatchesInSuperPatch)) {
      getline(fin,line);
      words = Common::StringOps::getWords(line);

      vector<std::string>::const_iterator itr = words.begin();
      for(;itr != words.end(); ++itr){
        _superPatch[i].getPatchIDs()[foundIDs] = StringOps::from_str<CFuint>(*itr);
        ++foundIDs;
        }
      }
    }
    else {
      i = i-1;
      getline(fin,line);
      }
  }

  (*fhandle)->close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::convertBack(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  throw Common::NotImplementedException (FromHere(),"Converting back from Gmsh not implemented.\n");
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::writeContinuousElements(ofstream& fout)
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

  fout << "!GEOM_POLYORDER " << _order << "\n";
  fout << "!SOL_POLYORDER "  << _order << "\n";

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(_elementType[k].getNbNodesPerCell(),
                                      _elementType[k].getOrderPerCell(),
                                      _dimension) << "\n";
  }

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
    const CFuint nbCellsPerType   = _elementType[k].getNbCellsPerType();
    const CFuint nbNodesPerCell   = _elementType[k].getNbNodesPerCell();
    const CFuint nbStatesPerCell  = nbNodesPerCell;
    const CFuint typeID           = _elementType[k].getTypeID();

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j) {
        const CFuint elemIdx = _mapNodeIdxPerElemTypeTable[typeID][j];
        fout << _elementType[k].getTableConnectivity()(i,elemIdx) << " " ;
      }
      for (CFuint j = 0; j < nbStatesPerCell; ++j) {
        const CFuint elemIdx = _mapNodeIdxPerElemTypeTable[typeID][j];
        fout << _elementType[k].getTableConnectivity()(i,elemIdx) << " " ;
      }
      fout << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << _isWithSolution << "\n";

  if (_isWithSolution) {
    fout.precision(12);
    fout << *_variables;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_NODE " << "\n";
  for (CFuint k = 0; k < _nbUpdatableNodes; ++k) {
    for (CFuint j = 0; j < _dimension; ++j)
      {
  fout.precision(12);
  fout << (*_coordinate)(k,j) << " ";
      }
    fout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
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

  // Reserve mem in map not possible: don't know how many times a key will
  // appear
  // m_nodeElement.reserve (_nbCells);

  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  fout << "!GEOM_POLYORDER " << (CFuint) _order << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER0 << "\n";
  // this CFPolyOrder::ORDER0 can only be set if here we know that we are dealing
  // with CellCenterFEM

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(_elementType[k].getNbNodesPerCell(),
                                      _elementType[k].getOrderPerCell(),
                                      _dimension) << "\n";
  }

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
    const CFuint typeID         = _elementType[k].getTypeID();

    // For use in the TRS
    _elementType[k].setCellStartIDX (countElem);

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j)
      {
        const CFuint elemIdx = _mapNodeIdxPerElemTypeTable[typeID][j];
        const CFuint nodeid = _elementType[k].getTableConnectivity()(i,elemIdx);
        fout << nodeid << " " ;
        m_nodeElement.insert(nodeid, countElem);
      }
      fout << countElem << "\n"; // cellID == stateID

      ++countElem;
    }

  }

  // sort map
  m_nodeElement.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
  fout << "!LIST_STATE " << _isWithSolution << "\n";

  CFuint nbVariables = PhysicalModelStack::getActive()->getNbEq();

  std::valarray<CFreal> averageState(nbVariables);

  /// @todo the following needs to be tested
  if (_isWithSolution) {
    for (CFuint k = 0; k < getNbElementTypes(); ++k) {
      const CFuint nbCellsPerType = _elementType[k].getNbCellsPerType();
      const CFuint nbNodesPerCell = _elementType[k].getNbNodesPerCell();

      cf_assert(nbCellsPerType > 0);
      for (CFuint i = 0; i < nbCellsPerType; ++i) {
         averageState = 0.0;
         for (CFuint j = 0; j < nbNodesPerCell; ++j) {
          const CFuint nodeID = _elementType[k].getTableConnectivity()(i,j);
          const vector<CFreal>& nodalState = _variables->getRow(nodeID);
          for(CFuint v = 0; v < nbVariables; ++v) {
            averageState[v] += nodalState[v];
          }
        }
        averageState /= nbNodesPerCell; // compute the average state
        for (CFuint iVar = 0; iVar < nbVariables; ++iVar) {
          fout.precision(12);
          fout << averageState[iVar] << " ";
        }
        fout << "\n";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////


void Gmsh2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;

  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < _nbPatches; ++iPatch) {
    const CFuint codeID = _patch[iPatch].getPatchCode();
    mapCP[codeID] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << _nbSuperPatches << " " << nbUnsetTRS << "\n";
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    const CFuint nbTRsInTRS = _superPatch[iTRS].getNbPatchesInSuperPatch();
    const std::string nameTRS  = _superPatch[iTRS].getSuperPatchName();

    fout << "!TRS_NAME " << nameTRS << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
      const CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
      const CFuint curPatch = mapCP.find(patchID)->second;
      fout << _patch[curPatch].getNbFacesInPatch() << " ";
    }
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

      const CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
      const CFuint curPatch = mapCP.find(patchID)->second;
     
      const CFuint nbFacesInPatch = _patch[curPatch].getNbFacesInPatch();
     

      for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
      {
    const CFuint nbNodesPerFace = _patch[curPatch].getFaceData()[iFace].getNbNodesInFace();
    const CFuint nbStatesPerFace = nbNodesPerFace;
    fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
    for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
      fout << _patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
      }
    for (CFuint iState = 0; iState < nbStatesPerFace; ++iState) {
      fout << _patch[curPatch].getFaceData()[iFace].getFaceNodes()[iState] << " ";
    }
    fout << "\n";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
{
   CFAUTOTRACE;

   typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;
   MapCodeNbPatch mapCP;

   for (CFuint iPatch = 0; iPatch < _nbPatches; ++iPatch) {
      const CFuint codeP = _patch[iPatch].getPatchCode();
      mapCP[codeP] = iPatch;
   }

   //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
   //fout << "!NB_TRSs " << _nbSuperPatches << " " << nbUnsetTRS << "\n";
   fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

   //only boundary TRS are listed !!!
   for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

      const CFuint nbTRsInTRS = _superPatch[iTRS].getNbPatchesInSuperPatch();
      const std::string nameTRS  = _superPatch[iTRS].getSuperPatchName();

      fout << "!TRS_NAME " << nameTRS << "\n";
      fout << "!NB_TRs "<< nbTRsInTRS << "\n";
      fout << "!NB_GEOM_ENTS ";

      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
         const CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
         const CFuint curPatch = mapCP.find(patchID)->second;
         fout << _patch[curPatch].getNbFacesInPatch() << " ";
      }
      fout << "\n";
      fout << "!GEOM_TYPE Face" << "\n";
      fout << "!LIST_GEOM_ENT" << "\n";

      // Vector holding element numbers for cell lookup
      vector<CFuint> elemens;
      /// Vector holding the nodes we care about
      vector<CFuint> trNodes;
      vector<CFuint> eleNodes;

      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR)
      {

         const CFuint patchID = _superPatch[iTRS].getPatchIDs()[iTR];
         const CFuint curPatch = mapCP.find(patchID)->second;
         const CFuint nbFacesInPatch = _patch[curPatch].getNbFacesInPatch();

         for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
         {
            const CFuint nbNodesPerFace = _patch[curPatch].getFaceData()[iFace].getNbNodesInFace();
            const CFuint nbStatesPerFace = 1;
            fout << nbNodesPerFace << " " << nbStatesPerFace << " ";

            elemens.clear();
            /// Find the corresponding cell ID
            typedef CFMultiMap<CFuint,CFuint>::MapIterator Iter;
            pair<Iter,Iter> res;

            trNodes.clear();

            for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
            {
               const CFuint nodeid = _patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode];
               /// @note the mapping between node local indexes in COOLFluiD and Gmsh is not used here,
               /// the face type is not immediately available. Should be okay, since the numbering of the nodes
               /// is the same for 1D and 2D elements in COOLFluiD and Gmsh
               fout << nodeid << " " ;

               trNodes.push_back(nodeid);

               // Store all the elements that contain this node in elemns
               bool fo = false;
	       res = m_nodeElement.find(nodeid, fo);
	       if (!fo) cout << "Gambit2CFmesh::writeDiscontinuousTrsData() => node " << nodeid  << " not found!\n";
	       for (Iter i = res.first; i != res.second; ++i) {
		 elemens.push_back(i->second);
               }
	       
            }
	    
            // Sort elemns
            sort(elemens.begin(), elemens.end());

            // sort nodes
            sort(trNodes.begin(), trNodes.end());

            // Find an element that contains all our nodes
            CFuint count = 0;
            CFuint last = elemens.front();
            for (CFuint i = 0; i < elemens.size() && count != nbNodesPerFace; ++i)
            {
               if (elemens[i] == last)
               {
                ++count;
                continue;
               }
               // New type
               count=1;
               last = elemens[i];
            }

            cf_assert (count <= nbNodesPerFace);

#ifdef NDEBUG
            // Double check the algorithm above
            vector<CFuint> uniq_ele (elemens);
            uniq_ele.erase(std::unique(uniq_ele.begin(), uniq_ele.end()),
                  uniq_ele.end());
            for (CFuint i=0; i<uniq_ele.size(); ++i)
            {
              CFuint mycount = std::count (elemens.begin(), elemens.end(), uniq_ele[i]);

               if (mycount >= nbNodesPerFace)
               {
                  if (count != nbNodesPerFace)
                  {
                     cerr << "For " << iTRS << "," << iTR << "," << iFace <<
                        ": count=" << count << ", mycount=" << mycount <<
                        ", nbNodesPerFace=" << nbNodesPerFace << endl;
                  }
               }
            }

#endif
            if (count != nbNodesPerFace)
               cerr << "No match found for " << iTRS << "," << iTR << ","
                  << iFace << " (count=" << count << ", nbNodesPerFace="
                  << nbNodesPerFace << ")" << endl;

            fout << last << "\n";

            // Verify that the element contains the state & nodes

            // Find element type
            CFuint eletype = getNbElementTypes();
            CFuint inCell = static_cast<CFuint>(-1);
            for (CFint e = getNbElementTypes()-1; e >= 0; --e)
            {
               if (last >= _elementType[e].getCellStartIDX ())
               {
                  eletype = e;
                  inCell = last - _elementType[e].getCellStartIDX();
                  break;
               }
            }

            cf_assert (eletype < getNbElementTypes());
            cf_assert (inCell < _elementType[eletype].getNbCellsPerType());

            const CFuint nbNodes =
               _elementType[eletype].getNbNodesPerCell();

            eleNodes.clear();
            for (CFuint i = 0; i < nbNodes; ++i)
               eleNodes.push_back
                  (_elementType[eletype].getTableConnectivity()(inCell,i));

            sort(eleNodes.begin(), eleNodes.end());

            const bool error = (!includes(eleNodes.begin(), eleNodes.end(),
                                    trNodes.begin(), trNodes.end()));

            if (error)
            {
               cerr << "Error in TR state matching:\n"
                  << "Selected TRS,TR,GEO: " << iTRS << "," << iTR << ","
                  << iFace << "\n";
               cerr << "Matched element: (and state)" << last <<
                  ", type=" << eletype << "\n";
               cerr << "Geo nodes: ";
               for (CFuint i = 0; i < trNodes.size(); ++i)
                  cerr << trNodes[i] << " ";
               cerr << endl;
               cerr << "Element nodes: ";
               for (CFuint i = 0; i<eleNodes.size(); ++i)
                  cerr << eleNodes[i] << " ";
               cerr << endl;
               cerr << "Possible elements: ";
               for (CFuint i=0; i<elemens.size(); ++i)
                  cerr << elemens[i] << " ";

               cerr << endl;
               cerr << endl;


            }

            fout.flush ();
         }
      }
   }
}

//////////////////////////////////////////////////////////////////////////////

void Gmsh2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  if(!_isFileRead) {
    try {
      checkFormat(filepath);
      if (_fileFormatVersion == 1)
      {
        readSPInnerCells(filepath);
        readGmshFileVersion1(filepath);
        readSPFile(filepath);
      }
      else if (_fileFormatVersion == 2)
      {
        readGmshFileVersion2(filepath);
      }
      _isFileRead = true;
    }
    catch (Common::Exception& e) {
      CFout << e.what() << "\n";
      throw;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Gmsh2CFmesh

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
