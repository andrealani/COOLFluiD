// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFMap.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Dpl2CFmesh/Dpl2CFmeshConverter.hh"
#include "Dpl2CFmesh/Dpl2CFmesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Dpl2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Dpl2CFmeshConverter,
               MeshFormatConverter,
               Dpl2CFmeshModule,
               1>
Dpl2CFmeshConverterProvider("Dpl2CFmesh");

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("isHybrid","Flag to tell if you convert an hybrid mesh");
}

//////////////////////////////////////////////////////////////////////////////

Dpl2CFmeshConverter::Dpl2CFmeshConverter (const std::string& name)
: MeshFormatConverter(name),
  _dimension(0),
  _nbCells(0),
  _nbUpdatableNodes(0),
  _nbFaces(0),
  _nbPatches(0),
  _nbSuperPatches(0),
  _elementType(0),
  _patch(0),
  _superPatch(0),
  _update(static_cast<CFuint>(0),static_cast<CFuint>(0)),
  _coordinate(CFNULL),
  _variables(CFNULL),
  _isWithSolution(false),
  _isFileRead(false),
  _nbUpdatableStates(_nbUpdatableNodes),
  _nodesPerElemTypeTable(2),
  _isHybrid(false)
{
   addConfigOptionsTo(this);
  /// Build the nbNodes per ElemTypeTable
  _nodesPerElemTypeTable[0] = 3;
  _nodesPerElemTypeTable[1] = 4;

  _isHybrid = false;
   setParameter("isHybrid",&_isHybrid);

}

//////////////////////////////////////////////////////////////////////////////

Dpl2CFmeshConverter::~Dpl2CFmeshConverter()
{
  deletePtr(_coordinate);
  deletePtr(_variables);
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::checkFormat(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path meshFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
  path meshFile = change_extension(filepath, getOriginExtension());
#endif

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(meshFile);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;

  // check
  getDplWordsFromLine(fin,line,lineNb,words);
  if(words.size() != 3) callDplFileError("Malformed file format.",lineNb,meshFile);
  if(words[0] != "unstructured") callDplFileError("Malformed file format",lineNb,meshFile);
  if(words[1] != "grid") callDplFileError("Malformed file format",lineNb,meshFile);
  if(words[2] != "data") callDplFileError("Malformed file format",lineNb,meshFile);

  // check nb Elements,
  // nb nodes, nb Boundary Faces
  // nb patches
  getDplWordsFromLine(fin,line,lineNb,words);
  if(words.size() != 3) callDplFileError("Wrong number of parameters.",lineNb,meshFile);
  const CFint nbElements = StringOps::from_str<CFint>(words[0]);
  if(nbElements < 0) callDplFileError("Negative number of Elements",lineNb,meshFile);
  const CFint nbNodes = StringOps::from_str<CFint>(words[1]);
  if(nbNodes < 0) callDplFileError("Negative number of Nodes",lineNb,meshFile);

  ///@todo add many more checks

  // Close file
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::callDplFileError(const std::string& msg,
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

void Dpl2CFmeshConverter::getDplWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::readDplFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path meshFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
   path meshFile = change_extension(filepath, getOriginExtension());
#endif

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(meshFile);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;

  // read general information
  // check
  getDplWordsFromLine(fin,line,lineNb,words);
  if(words.size() != 3) callDplFileError("Malformed file format.",lineNb,meshFile);
  if(words[0] != "unstructured") callDplFileError("Malformed file format",lineNb,meshFile);
  if(words[1] != "grid") callDplFileError("Malformed file format",lineNb,meshFile);
  if(words[2] != "data") callDplFileError("Malformed file format",lineNb,meshFile);
//
  // Read Number of nodes
  getDplWordsFromLine(fin,line,lineNb,words);
  _nbUpdatableNodes = StringOps::from_str<CFint>(words[1]);

  // Read Number of Elements
  _nbCells = StringOps::from_str<CFint>(words[0]);

  ///Dpl is always 2D.
  _dimension = 2;

  CFuint nbElementTypes;
  if(_isHybrid) nbElementTypes = 2;
  else nbElementTypes = 1;

  // Elements are not grouped by element types...
  // First loop over all the elements to know the number of element types
  // Also get the number of Physical Regions (<->SP) and elemnt Regions (<->Patches)

  // Read and allocate element type information
  _elementType.resize(nbElementTypes);

  // Build some temporaries
  CFuint nbNodes;

  // Temporary int to count the number of elements
  CFuint totalNbElements;
  totalNbElements = 0;
  std::vector<CFuint> nbElemPerType(nbElementTypes);

  for (CFuint i = 0; i < _nbCells; ++i){
    getDplWordsFromLine(fin,line,lineNb,words);
    nbNodes = StringOps::from_str<CFint>(words[0]);

    nbElemPerType[nbNodes-3] += 1;
    totalNbElements += 1;
  }

  // Now we know the number of:
  // - Elements
  // - ElementTypes and Nb of Elements per Type
  /// This allows the allocation of the memory

  // For each elementType present in the mesh, set the typeID
  for (CFuint iType = 0; iType < nbElementTypes; ++iType){
    _elementType[iType].setTypeID(iType);
  }

  // For each elementType present in the mesh, set the nbNodesPerCell and NbCellsPerType
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType) {
    _elementType[iType].setNbNodesPerCell(_nodesPerElemTypeTable[_elementType[iType].getTypeID()]);
    _elementType[iType].setNbCellsPerType(nbElemPerType[_elementType[iType].getTypeID()]);
    _elementType[iType].setCurrentIndex(0);
  }

  // Create storage space for the connectivity tables
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    _elementType[iType].createCellNodeConnectivity();
  }

  /// Re-loop over the elements to effectively read the elements
  //Change the line Number to reloop over the elements
  lineNb -= _nbCells;

  /// Go back to the beginning of the elements
  /// @todo try to do this without closing/opening...

  fhandle->close();
  fhandle->open(meshFile);

  for (CFuint i = 0; i < lineNb; ++i){
    getline(fin,line);
    }

  CFuint elemTypeIdx = 0;
  CFuint inTypeIdx;

  /// Reading the Elements
  for (CFuint i = 0; i < _nbCells; ++i){
    getDplWordsFromLine(fin,line,lineNb,words);

    //What is the type of the cell??
    nbNodes = StringOps::from_str<CFint>(words[0]);

    std::vector<CFuint> nodesIdx(nbNodes);
    for(CFuint iNode=0;iNode < nbNodes; iNode++)
    {
      nodesIdx[iNode] = StringOps::from_str<CFint>(words[iNode+1]);
    }

    elemTypeIdx = nbNodes - 3;

    inTypeIdx = _elementType[elemTypeIdx].getCurrentIndex();

    for (CFuint j = 0; j < _elementType[elemTypeIdx].getNbNodesPerCell(); ++j){
      _elementType[elemTypeIdx].getTableConnectivity()(inTypeIdx,j) = nodesIdx[j] - 1;
    }
    _elementType[elemTypeIdx].setCurrentIndex(inTypeIdx+1);

  }


  /// Reading the Coordinates
  //pass two lines
  getDplWordsFromLine(fin,line,lineNb,words);
  getDplWordsFromLine(fin,line,lineNb,words);

  // read node coordinates
  _coordinate = new Table<CFreal>(_nbUpdatableNodes, _dimension);
  for (CFuint iNode = 0; iNode < _nbUpdatableNodes; ++iNode)
  {
    getDplWordsFromLine(fin,line,lineNb,words);
    for (CFuint iDim = 0; iDim < _dimension; ++iDim){
      (*_coordinate)(iNode,iDim) = StringOps::from_str<CFreal>(words[iDim]);
      }
  }

  /// Reading the Boundaries/Patches
  // Read the number of boundaries/patches
  getDplWordsFromLine(fin,line,lineNb,words);

  _nbPatches = StringOps::from_str<CFuint>(words[0]);
  // allocate patch storage space
  _patch.resize(_nbPatches);

  for(CFuint iBoundary=0; iBoundary < _nbPatches ; iBoundary++)
  {
    getDplWordsFromLine(fin,line,lineNb,words);
    CFuint nbFacesInPatch = StringOps::from_str<CFuint>(words[0]);
    CFuint patchCode = StringOps::from_str<CFuint>(words[1]);
    _patch[iBoundary].setPatchCode(patchCode);
    _patch[iBoundary].setNbFacesInPatch(nbFacesInPatch);

    for(CFuint iElem=0; iElem < nbFacesInPatch ; iElem++)
    {
      getDplWordsFromLine(fin,line,lineNb,words);
      CFuint  nbFaceNodes = 2;
// std::cout << "    CellID: " << words[nbFaceNodes*2].StringOps::from_str<CFuint>() - 1 << std::endl;
// std::cout << "    Words: " << line << std::endl;
      _patch[iBoundary].getFaceData()[iElem].setCellID(StringOps::from_str<CFuint>(words[nbFaceNodes*2]) - 1);
      _patch[iBoundary].getFaceData()[iElem].setNbNodesInFace(nbFaceNodes);
      for (CFuint j = 0; j < nbFaceNodes; ++j) {
        _patch[iBoundary].getFaceData()[iElem].getFaceNodes()[j] = StringOps::from_str<CFuint>(words[j]) - 1;
      }
    }
  }

  fhandle->close();


}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::readSPFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path fileSP = boost::filesystem::path(filepath).replace_extension(".SP");
#else
  path fileSP = change_extension(filepath,".SP");
#endif
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(fileSP);

  std::string line;
  vector<std::string> words;

  getline(fin,line);
  words = Common::StringOps::getWords(line);

  if(words.size() != 1) {
    throw BadFormatException (FromHere(),"Bad number of parameters in file " + fileSP.string());
  }

  CFuint nbSuperPatches = StringOps::from_str<CFuint>(words[0]);

  if(nbSuperPatches > _nbPatches) {
    throw BadFormatException (FromHere(),"Number of SuperPatches bigger than Patches " + fileSP.string());
  }

  _superPatch.resize(nbSuperPatches);

  for (CFuint i = 0; i < nbSuperPatches; ++i) {

    getline(fin,line);
    words = Common::StringOps::getWords(line);

    if(words.size() != 2) {
      throw BadFormatException (FromHere(),"Bad number of parameters in file " + fileSP.string());
    }

    std::string superPatchName = words[0];
    CFint nbPatchesInSuperPatch = StringOps::from_str<CFint>(words[1]);

    if((CFuint)nbPatchesInSuperPatch > _nbPatches ||
       nbPatchesInSuperPatch <= 0) {
      throw BadFormatException (FromHere(),"Bad number of Patches: "
             + words[1]
             + " in SuperPatch "
             + fileSP.string());
    }

    _superPatch[i].setSuperPatchName(superPatchName);
    _superPatch[i].setNbPatchesInSuperPatch(nbPatchesInSuperPatch);

    CFuint foundIDs = 0;
    while (foundIDs < (CFuint)nbPatchesInSuperPatch) {
      getline(fin,line);
      words = Common::StringOps::getWords(line);

      vector<std::string>::const_iterator itr = words.begin();
      for(;itr != words.end(); ++itr){
  _superPatch[i].getPatchIDs()[foundIDs] = StringOps::from_str<CFuint>(*itr);
  ++foundIDs;
      }
    }
  }

  fhandle->close();
}


//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::writeContinuousElements(ofstream& fout)
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
    CFuint nbCellsPerType = _elementType[k].getNbCellsPerType();
    CFuint nbNodesPerCell = _elementType[k].getNbNodesPerCell();
    CFuint nbStatesPerCell  = nbNodesPerCell;

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j) {
        fout << _elementType[k].getTableConnectivity()(i,j) <<" " ;
      }
      for (CFuint j = 0; j < nbStatesPerCell; ++j) {
        fout << _elementType[k].getTableConnectivity()(i,j) <<" " ;
      }
      fout << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << _isWithSolution << "\n";

  if (_isWithSolution) {
    fout.precision(12);
    fout << *_variables;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::writeNodes(ofstream& fout)
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

void Dpl2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
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

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j) {
        fout << _elementType[k].getTableConnectivity()(i,j) << " " ;
      }
      fout << countElem << "\n"; // cellID == stateID
      ++countElem;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
  fout << "!LIST_STATE " << _isWithSolution << "\n";

  CFuint nbVariables = PhysicalModelStack::getActive()->getNbEq();

  std::valarray<CFreal> averageState(nbVariables);

  if (_isWithSolution) {
    for (CFuint k = 0; k < getNbElementTypes(); ++k) {
      CFuint nbCellsPerType = _elementType[k].getNbCellsPerType();
      CFuint nbNodesPerCell = _elementType[k].getNbNodesPerCell();

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

void Dpl2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;

  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < _nbPatches; ++iPatch) {
    CFuint codeID = _patch[iPatch].getPatchCode();
    mapCP[codeID] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << _nbSuperPatches << " " << nbUnsetTRS << "\n";
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = _superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = _superPatch[iTRS].getSuperPatchName();

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

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
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
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = _superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = _superPatch[iTRS].getSuperPatchName();

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

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  if(!_isFileRead) {
    try {
      readDplFile(filepath);
      readSPFile(filepath);
      _isFileRead = true;
    }
    catch (Common::Exception& e) {
      CFout << e.what() << "\n";
      throw;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Dpl2CFmeshConverter::convertBack(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  throw Common::NotImplementedException (FromHere(),"Converting back from Dpl not implemented.\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Dpl2CFmesh

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
