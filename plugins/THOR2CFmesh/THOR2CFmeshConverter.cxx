// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>

#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Common/Stopwatch.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "THOR2CFmesh/THOR2CFmeshConverter.hh"
#include "THOR2CFmesh/CheckNodeNumberingTetra.hh"
#include "THOR2CFmesh/CheckNodeNumberingHexa.hh"
#include "THOR2CFmesh/CheckNodeNumberingPyram.hh"
#include "THOR2CFmesh/CheckNodeNumberingPrism.hh"
#include "THOR2CFmesh/THOR2CFmesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<THOR2CFmeshConverter,
               MeshFormatConverter,
               THOR2CFmeshModule,
               1>
thor2CFmeshConverterProvider("THOR2CFmesh");

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SolutionOrder","Order of the solution space to be created in the converted file.");
}

//////////////////////////////////////////////////////////////////////////////

THOR2CFmeshConverter::THOR2CFmeshConverter (const std::string& name)
: MeshFormatConverter(name),
  m_dimension(0),
  m_nbLamVariables(0),
  m_nbTurbVariables(0),
  m_nbCells(0),
  m_nbUpdatableNodes(0),
  m_nbFaces(0),
  m_nbPatches(0),
  m_elementType(0),
  m_patch(0),
  m_superPatch(0),
  m_update(static_cast<CFuint>(0),static_cast<CFuint>(0)),
  m_coordinate(CFNULL),
  m_variables(CFNULL),
  m_scon(CFNULL),
  m_isWithSolution(false),
  m_isFileRead(false),
  m_offsetTHOR(true),
  m_nbUpdatableStates(m_nbUpdatableNodes)
{
  addConfigOptionsTo(this);

  m_solOrderStr = "P1";
  setParameter( "SolutionOrder",      &m_solOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

THOR2CFmeshConverter::~THOR2CFmeshConverter()
{
  deletePtr(m_coordinate);
  deletePtr(m_variables);
  deletePtr(m_scon);
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::configure ( Config::ConfigArgs& args )
{
  MeshFormatConverter::configure(args);

  m_solOrder = CFPolyOrder::Convert::to_enum( m_solOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::checkFormat(const boost::filesystem::path& filepath)
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

  // check dimension, nb laminar variables, nb turb variables
  getTHORWordsFromLine(fin,line,lineNb,words);
  if(words.size() != 3) callTHORFileError("Wrong number of parameters.",lineNb,meshFile);
  const CFint dim     = Common::StringOps::from_str<CFint>(words[0]);
  const CFint lamVars = Common::StringOps::from_str<CFint>(words[1]);
  const CFint trbVars = Common::StringOps::from_str<CFint>(words[2]);
  if(dim != 2 && dim != 3) callTHORFileError("Malformed dimension",lineNb,meshFile);
  if(lamVars < 0) callTHORFileError("Wrong number of laminar variables",lineNb,meshFile);
  if(trbVars < 0) callTHORFileError("Wrong number of turbulent variables",lineNb,meshFile);

  // check nb Elements,
  // nb nodes, nb Boundary Faces
  // nb patches
  getTHORWordsFromLine(fin,line,lineNb,words);

  if(words.size() != 4) callTHORFileError("Wrong number of parameters.",lineNb,meshFile);
  const CFint nbElems = Common::StringOps::from_str<CFint>(words[0]);
  const CFint nbNodes = Common::StringOps::from_str<CFint>(words[1]);
  const CFint nbFaces = Common::StringOps::from_str<CFint>(words[2]);
  const CFint nbPatch = Common::StringOps::from_str<CFint>(words[3]);
  if(nbElems < 0) callTHORFileError("Negative number of Elements",lineNb,meshFile);
  if(nbNodes < 0) callTHORFileError("Negative number of Nodes",lineNb,meshFile);
  if(nbFaces < 0) callTHORFileError("Negative number of Faces",lineNb,meshFile);
  if(nbPatch < 0) callTHORFileError("Negative number of Patches",lineNb,meshFile);

  // this arrays store a check if the node/states are used by some elements
  // to avoid having loose nodes on the mesh
  std::valarray<bool> nodeUsed(nbNodes);
  nodeUsed = false;

  // check nb Elements Types,
  getTHORWordsFromLine(fin,line,lineNb,words);

  if(words.size() != 1) callTHORFileError("Wrong number of parameters.",lineNb,meshFile);
  const CFint nbElementTypes = Common::StringOps::from_str<CFint>(words[0]);
  if(nbElementTypes < 0) callTHORFileError("Negative number of Element Types",lineNb,meshFile);
  if(dim == DIM_3D && nbElementTypes > 4)
    callTHORFileError("Number of Element Types greater than 4 in 3D mesh",lineNb,meshFile);
  if(dim == DIM_2D && nbElementTypes > 2)
    callTHORFileError("Number of Element Types greater than 2 in 2D mesh",lineNb,meshFile);

  // check number cells per type
  CFuint totalNbElems = 0;
  Table<CFint> elemTypes(nbElementTypes,2);
  for (CFint k = 0; k < nbElementTypes; ++k) {
    getTHORWordsFromLine(fin,line,lineNb,words);

    if (words.size() != 2)
      callTHORFileError("Wrong number of parameters.",lineNb,meshFile);
    const CFint nbNodesPerCell = Common::StringOps::from_str<CFint>(words[0]);
    const CFint nbCellsInType  = Common::StringOps::from_str<CFint>(words[1]);
    if(dim == DIM_2D && nbNodesPerCell != 3 && nbNodesPerCell != 4)
       callTHORFileError("Wrong number of nodes per cell in 2D",lineNb,meshFile);
    if(dim == DIM_3D && nbNodesPerCell != 4 &&
       nbNodesPerCell != 5 &&
       nbNodesPerCell != 6 &&
       nbNodesPerCell != 8)
      callTHORFileError("Wrong number of nodes per cell in 3D",lineNb,meshFile);

    elemTypes(k,0) = nbNodesPerCell;
    elemTypes(k,1) = nbCellsInType;
    totalNbElems += nbCellsInType;
  }

  // check that the sum of elements per type equals the total nb elements
  if(totalNbElems != (CFuint)nbElems) {
    callTHORFileError("Bad sum of number of cells per type",lineNb,meshFile);
  }

  // check per element type connectivity
  for (CFint k = 0; k < nbElementTypes; ++k) {
    const CFuint nbNodesPerCell = elemTypes(k,0);
    const CFuint nbCellsInType = elemTypes(k,1);
    for (CFuint i = 0; i < nbCellsInType; ++i) {
      getTHORWordsFromLine(fin,line,lineNb,words);

      if(words.size() != nbNodesPerCell)
  callTHORFileError("Wrong number of nodes in cell",lineNb,meshFile);
      for (CFuint n = 0; n < nbNodesPerCell; ++n) {
  const CFint idx = Common::StringOps::from_str<CFint>(words[n]);
  if(idx < 1 || idx > nbNodes) {
    callTHORFileError("Bad nodex index",lineNb,meshFile);
  }
  nodeUsed[idx-1] = true;
      }
    }
  }

  // check patches
  for (CFint k = 0; k < nbPatch; ++k) {
    getTHORWordsFromLine(fin,line,lineNb,words);
    if(words.size() != 2)
      callTHORFileError("Wrong number of parameters",lineNb,meshFile);

    const CFint nbFacesInPatch = Common::StringOps::from_str<CFint>(words[1]);
    if(nbFacesInPatch < 0 || nbFacesInPatch > nbFaces)
      callTHORFileError("Wrong number of parameters",lineNb,meshFile);
    for (CFint f = 0; f < nbFacesInPatch; ++f) {
      getTHORWordsFromLine(fin,line,lineNb,words);
      if(words.size() <= 2 || words.size() > 6)
        callTHORFileError("Wrong number of parameters",lineNb,meshFile);

      const CFint elemIdx = Common::StringOps::from_str<CFint>(words[0]);
      if(elemIdx < 1 || elemIdx > nbElems)
        callTHORFileError("Bad element index",lineNb,meshFile);

      const CFint nbNodesInFace = Common::StringOps::from_str<CFint>(words[1]);

      if(dim == 2 && nbNodesInFace != 2)
        callTHORFileError("Wrong number of faces in 2D",lineNb,meshFile);

      if(dim == 3 && nbNodesInFace != 3 && nbNodesInFace != 4)
        callTHORFileError("Wrong number of faces in 3D",lineNb,meshFile);

      for (CFint n = 0; n < nbNodesInFace; ++n) {
        const CFint idx = Common::StringOps::from_str<CFint>(words[n+2]);
          if(idx < 1 || idx > nbNodes)
            callTHORFileError("Bad nodex index",lineNb,meshFile);
      }
    }
  }

  /// @todo More checks to THOR Format files need to be added.
  ///       Check Node listing
  ///       Check Coordinates and solution

  // check that all nodes are used
  bool nodesLoose = false;
  for (CFuint i = 0; i < (CFuint)nbNodes; ++i) {
    if(nodeUsed[i] != true) {
      nodesLoose = true;
      CFerr << "Loose Node: " << i << "\n";
    }
  }
  if(nodesLoose) {
    throw BadFormatException (FromHere(),"Loose nodes exist in THOR file");
  }

  fhandle->close();


}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::callTHORFileError(const std::string& msg,
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

void THOR2CFmeshConverter::getTHORWordsFromLine(ifstream& fin,
            std::string& line,
            CFuint&  lineNb,
            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::readTHOR(const boost::filesystem::path& filepath)
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

  CFuint nbElementTypes;

  // read general information
  fin >> m_dimension >> m_nbLamVariables >> m_nbTurbVariables;

  fin >> m_nbCells >> m_nbUpdatableNodes >> m_nbFaces >> m_nbPatches;
  fin >> nbElementTypes;

  //read and allocate element type information
  m_elementType.resize(nbElementTypes);

  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    CFuint  nbNodesPerCell;
    CFuint nbCellsPerType;
    fin >> nbNodesPerCell >> nbCellsPerType;
    m_elementType[k].setNbNodesPerCell(nbNodesPerCell);
    m_elementType[k].setNbCellsPerType(nbCellsPerType);
  }

  // create storage space for the connectivity tables
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    m_elementType[k].createCellNodeConnectivity();
  }

  // read per element type connectivity
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    for (CFuint i = 0; i <  m_elementType[k].getNbCellsPerType(); ++i) {
      for (CFuint j = 0; j < m_elementType[k].getNbNodesPerCell(); ++j) {
        fin >> m_elementType[k].getTableConnectivity()(i,j);
      }
    }
  }

  // allocate patch storage space
  m_patch.resize(m_nbPatches);

  // read and allocate patch face information
  for (CFuint k = 0; k < m_nbPatches; ++k) {
    CFuint patchCode;
    CFuint nbFacesInPatch;
    fin >> patchCode >> nbFacesInPatch;

    m_patch[k].setPatchCode(patchCode);
    m_patch[k].setNbFacesInPatch(nbFacesInPatch);

    for (CFuint i = 0; i < nbFacesInPatch; ++i) {
      CFuint  nbFaceNodes;
      CFuint  cellID;
      fin >> cellID >> nbFaceNodes;
      m_patch[k].getFaceData()[i].setCellID(cellID);
      m_patch[k].getFaceData()[i].setNbNodesInFace(nbFaceNodes);
      for (CFuint j = 0; j < nbFaceNodes; ++j) {
        fin >> m_patch[k].getFaceData()[i].getFaceNodes()[j];
      }
    }
  }

  // read list of updatable nodes
  m_update.resize(m_nbUpdatableNodes);
  for (CFuint i = 0; i < m_nbUpdatableNodes; ++i)
    fin >> m_update[i];

  // read node coordinates
  m_coordinate = new Table<CFreal>(m_nbUpdatableNodes, m_dimension);
  for (CFuint i = 0; i < m_dimension; ++i)
    for (CFuint j = 0; j < m_nbUpdatableNodes; ++j)
      fin >> (*m_coordinate)(j,i);

  // allocate and read nodal solution
  CFuint nbVariables = m_nbLamVariables + m_nbTurbVariables;
  if (nbVariables > 0) {
    m_isWithSolution = true;

    m_variables = new Table<CFreal>(m_nbUpdatableNodes, nbVariables);

    for (CFuint i = 0; i < m_nbLamVariables; ++i)
      for (CFuint j = 0; j < m_nbUpdatableNodes; ++j)
  fin >> (*m_variables)(j,i);

    for (CFuint i = m_nbLamVariables; i < m_nbLamVariables + m_nbTurbVariables; ++i)
      for (CFuint j = 0; j < m_nbUpdatableNodes; ++j)
  fin >> (*m_variables)(j,i);
  }


// this is useful to see what patches
// uncomment for debuging purposes
#if 0 // put 1 to insert into code
  std::string patchFilename = Environment::DirPaths::getInstance().getWorkingDir()
    + "Patch.out"
  ofstream fout(patchFilename.c_str());
  for (CFuint k = 0; k < m_nbPatches; ++k) {
    fout << "iPatch " << k+1 << "\n";
    CFuint nbPatchFaces = m_patch[k].nbFaces;
    for (CFuint i = 0; i < nbPatchFaces; ++i) {
      CFuint nbFaceNodes = m_patch[k].face[i].nbNodes;
      for (CFuint j = 0; j < nbFaceNodes; ++j) {
  CFuint nodeid = m_patch[k].face[i].node[j];
  for (CFuint iDim = 0; iDim < m_dimension; ++iDim) {
    fout << (*m_coordinate)(nodeid-1,iDim) << " ";
  }
  fout << ", ";
      }
      fout << "\n";
    }
  }
  fout.close();
#endif

  fhandle->close();


}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::readSP(const boost::filesystem::path& filepath)
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

  CFuint nbSuperPatches = Common::StringOps::from_str<CFuint>(words[0]);

  if(nbSuperPatches > m_nbPatches) {
    throw BadFormatException (FromHere(),"Number of SuperPatches bigger than Patches " + fileSP.string());
  }

  m_superPatch.resize(nbSuperPatches);

  for (CFuint i = 0; i < nbSuperPatches; ++i) {

    getline(fin,line);
    words = Common::StringOps::getWords(line);

    if(words.size() != 2) {
      throw BadFormatException (FromHere(),"Bad number of parameters in file " + fileSP.string());
    }

    std::string superPatchName = words[0];
    CFint nbPatchesInSuperPatch = Common::StringOps::from_str<CFint>(words[1]);

    if((CFuint)nbPatchesInSuperPatch > m_nbPatches ||
       nbPatchesInSuperPatch <= 0) {
      throw BadFormatException (FromHere(),"Bad number of Patches: "
             + words[1]
             + " in SuperPatch "
             + fileSP.string());
    }

    m_superPatch[i].setSuperPatchName(superPatchName);
    m_superPatch[i].setNbPatchesInSuperPatch(nbPatchesInSuperPatch);

    CFuint foundIDs = 0;
    while (foundIDs < (CFuint)nbPatchesInSuperPatch) {
      getline(fin,line);
      words = Common::StringOps::getWords(line);

      vector<std::string>::const_iterator itr = words.begin();
      for(;itr != words.end(); ++itr){
        m_superPatch[i].getPatchIDs()[foundIDs] = Common::StringOps::from_str<CFuint>(*itr);
        ++foundIDs;
      }
    }
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::convertBack(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  // only reads if not yet been read
  readFiles(filepath);
  adjustToTHORNodeNumbering();

  writeTHOR(filepath);
  writeSP(filepath);
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeSP(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path outFile = boost::filesystem::path(filepath).replace_extension(".SP");
#else
  path outFile = change_extension(filepath, ".SP");
#endif

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(outFile);

  fout << m_superPatch.size() << "\n";
  for (CFuint i = 0; i < m_superPatch.size(); ++i)
  {
    fout << m_superPatch[i].getSuperPatchName() << " " << m_superPatch[i].getNbPatchesInSuperPatch() << "\n";
    for (CFuint j = 0; j < m_superPatch[i].getNbPatchesInSuperPatch(); ++j)
      fout << m_superPatch[i].getPatchIDs()[j] << "\n";
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeTHOR(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path outFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
  path outFile = change_extension(filepath, getOriginExtension());
#endif

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(outFile);

  // write general mesh information
  fout << m_dimension << " "
       << m_nbLamVariables  << " "
       << m_nbTurbVariables << "\n";
  fout << m_nbCells << " "
       << m_nbUpdatableNodes << " "
       <<  m_nbFaces << " "
       <<  m_nbPatches << "\n";
  fout << getNbElementTypes() << "\n";

  // write element typr information
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
    fout << m_elementType[k].getNbNodesPerCell() << " "
         << m_elementType[k].getNbCellsPerType() << "\n";

  // write element type
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
    for (CFuint i = 0; i <  m_elementType[k].getNbCellsPerType(); ++i) {
      for (CFuint j = 0; j < m_elementType[k].getNbNodesPerCell(); ++j)
        fout << m_elementType[k].getTableConnectivity()(i,j) <<" " ;
      fout << "\n";
    }

  // write patch face information
  for (CFuint k = 0; k < m_nbPatches; ++k) {
    fout << m_patch[k].getPatchCode() << " " << m_patch[k].getNbFacesInPatch() << "\n";

    for (CFuint i = 0; i < m_patch[k].getNbFacesInPatch(); ++i) {
      fout << m_patch[k].getFaceData()[i].getCellID() << " "
     << m_patch[k].getFaceData()[i].getNbNodesInFace() << " ";

      for (CFuint j = 0; j < m_patch[k].getFaceData()[i].getNbNodesInFace(); ++j)
        fout << m_patch[k].getFaceData()[i].getFaceNodes()[j] << " ";
      fout << "\n";
    }
  }

  // write list of updatable nodes
  CFuint count = 0;
  for (CFuint iNode = 0; iNode < m_nbUpdatableNodes; ++iNode) {
    ++count;
    fout << m_update[iNode];
    // put 10 indexes in each line
    if(count == 10 || iNode == m_nbUpdatableNodes-1) {
      fout << "\n";
      count = 0;
    }
    else fout << " ";
  }

  // write nodal coordinates
  count = 0;
  for (CFuint i = 0; i < m_dimension; ++i) {
    for (CFuint j = 0; j < m_nbUpdatableNodes; ++j) {
      ++count;
      fout << (*m_coordinate)(j,i);
      // put 4 coordinates in each line
      if(count == 4 || (i == m_dimension-1 && j == m_nbUpdatableNodes-1)) {
        fout << "\n";
        count = 0;
      }
      else fout << " ";
    }
  }


  // write nodal solution
  CFuint nbVariables = m_nbLamVariables + m_nbTurbVariables;

  for (CFuint i = 0; i < m_nbLamVariables; ++i) {
    for (CFuint j = 0; j < m_nbUpdatableNodes; ++j) {
      ++count;
      fout << (*m_variables)(j,i);
      if(count == 4 || (i == m_nbLamVariables-1 && j == m_nbUpdatableNodes-1)) {
  fout << "\n";
  count = 0;
      }
      else fout << " ";
    }
  }

  for (CFuint i = m_nbLamVariables; i < nbVariables; ++i) {
    for (CFuint j = 0; j < m_nbUpdatableNodes; ++j) {
      ++count;
      fout << (*m_variables)(j,i);
      if(count == 4 || (i == nbVariables-1 && j == m_nbUpdatableNodes-1)) {
  fout << "\n";
  count = 0;
      }
      else fout << " ";
    }
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;

  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < m_nbPatches; ++iPatch) {
    CFuint codeID = m_patch[iPatch].getPatchCode();
    mapCP[codeID] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << m_nbSuperPatches << " " << nbUnsetTRS << "\n";
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = m_superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = m_superPatch[iTRS].getSuperPatchName();

    fout << "!TRS_NAME " << nameTRS << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      fout << m_patch[curPatch].getNbFacesInPatch() << " ";
    }
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      CFuint nbFacesInPatch = m_patch[curPatch].getNbFacesInPatch();

      for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
      {

        CFuint nbNodesPerFace = m_patch[curPatch].getFaceData()[iFace].getNbNodesInFace();
        CFuint nbStatesPerFace = nbNodesPerFace;

        fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
        for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
          fout << m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
        }
        for (CFuint iState = 0; iState < nbStatesPerFace; ++iState) {
          fout << m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iState] << " " ;
        }
        fout << "\n";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;
  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < m_nbPatches; ++iPatch) {
    CFuint codeP = m_patch[iPatch].getPatchCode();
    mapCP[codeP] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << m_nbSuperPatches << " " << nbUnsetTRS << "\n";
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = m_superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = m_superPatch[iTRS].getSuperPatchName();

    fout << "!TRS_NAME " << nameTRS << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      fout << m_patch[curPatch].getNbFacesInPatch() << " ";
    }
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      CFuint nbFacesInPatch = m_patch[curPatch].getNbFacesInPatch();

      for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
      {
        CFuint nbNodesPerFace = m_patch[curPatch].getFaceData()[iFace].getNbNodesInFace();

        CFuint nbStatesPerFace;
        if (m_solOrder == CFPolyOrder::ORDER0)
          { nbStatesPerFace = 1; }
        else // ORDER == CFPolyOrder::ORDER1
          {
          if (m_dimension == 2)
            {nbStatesPerFace = 2; }
          else
            {
            nbStatesPerFace = 3;
            }
          }

        fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
        for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
          fout << m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
        }

        if (m_solOrder == CFPolyOrder::ORDER0)
        {
          // StateID hardlinked with the CellID
          fout << m_patch[curPatch].getFaceData()[iFace].getCellID();
        }
        else // ORDER == CFPolyOrder::ORDER1
        {
          cf_assert(m_scon != CFNULL);
          Common::Table< std::pair<CFuint,CFuint> >& scon = *m_scon;

          CFuint cellID = m_patch[curPatch].getFaceData()[iFace].getCellID();
          for (CFuint lnodeID = 0; lnodeID < scon.nbCols(cellID); ++lnodeID)
          {
            // find which n
            for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
            {
              if (scon(cellID,lnodeID).first == m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode])
                fout << scon(cellID,lnodeID).second << " " ;
            }
          }
        }
        fout << "\n" ;

      } // end loop faces
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::adjustToCFmeshNodeNumbering()
{
  CFAUTOTRACE;

  // node numbering in THOR-file start from 1,
  // while in COOLFluiD format they start from 0
  if (m_offsetTHOR) {
    CFLogDebugMin(
    "Adjusting index numbering to CFmesh numbering: offset -1" << "\n");
    offsetNumbering(-1);
    m_offsetTHOR = false;
  }

  //if the mesh is hybrid you can have mismatches if you change the node
  //ordering of a face => to take into account
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType) {
    checkNodeNumberings(m_dimension,
			m_elementType[iType],
			m_coordinate);
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::offsetNumbering(const CFint& offset)
{
  // offset node indexes
  for (CFuint k = 0; k < m_nbUpdatableNodes; ++k) {
    m_update[k] += offset;
  }

  // offset node indexes in table of connectivity
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    m_elementType[k].getTableConnectivity() += offset;
  }

  // offset node indexes in Patch information
  for (CFuint k = 0; k < m_nbPatches; ++k) {

    CFuint nbPatchFaces = m_patch[k].getNbFacesInPatch();
    for (CFuint i = 0; i < nbPatchFaces; ++i) {
      m_patch[k].getFaceData()[i].getCellID() += offset;

      CFuint nbFaceNodes = m_patch[k].getFaceData()[i].getNbNodesInFace();
      for (CFuint j = 0; j < nbFaceNodes; ++j) {
        m_patch[k].getFaceData()[i].getFaceNodes()[j] += offset;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::adjustToTHORNodeNumbering()
{
  // node numbering in THOR-file start from 1,
  // while in COOLFluiD format they start from 0
  if (!m_offsetTHOR) {
    CFLogDebugMin(
    "Adjusting index numbering to THOR numbering: offset +1" << "\n");
    offsetNumbering(+1);
    m_offsetTHOR = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void
THOR2CFmeshConverter::checkNodeNumberings(const CFuint& dim,
            ElementTypeTHOR& elemType,
            const Table<CFreal> *const nodes) const
{
  auto_ptr<CheckNodeNumbering> checkNN;

  switch(dim) {
  case DIM_2D:
    switch(elemType.getNbNodesPerCell()) {
    case 3:
      break;
    case 4:
      break;
    default:
      std::string msg = std::string("Wrong number of nodes in 2D thor element: ") +
                     StringOps::to_str(elemType.getNbNodesPerCell());
      throw BadValueException (FromHere(),msg);
    }
    break;
  case DIM_3D:
    switch(elemType.getNbNodesPerCell()) {
    case 4:
      checkNN.reset(new CheckNodeNumberingTetra(nodes));
      break;
    case 5:
      checkNN.reset(new CheckNodeNumberingPyram(nodes));
      break;
    case 6:
      checkNN.reset(new CheckNodeNumberingPrism(nodes));
      break;
    case 8:
      checkNN.reset(new CheckNodeNumberingHexa(nodes));
      break;
    default:
      std::string msg = std::string("Wrong number of nodes in 3D thor element: ") +
                     StringOps::to_str(elemType.getNbNodesPerCell());
      throw BadValueException (FromHere(),msg);
    }
    break;
  default:
    std::string msg = std::string("Wrong dimension in thor file. Can only be 2D or 3D.") +
                   StringOps::to_str(elemType.getNbNodesPerCell());
    throw BadValueException (FromHere(),msg);
  }

  if (checkNN.get() != CFNULL) {
    checkNN->checkElementNodalNumbering(elemType.getTableConnectivity());
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  if(!m_isFileRead)
  {
    try
    {
      readTHOR(filepath);
      readSP(filepath);
      m_isFileRead = true;
    }
    catch (Common::Exception& e)
    {
      CFout << e.what() << "\n";
      throw;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeContinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "THOR2CFmeshConverter::writeContinuousElements() => start\n");
  
  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  fout << "!NB_NODES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  fout << "!NB_STATES "
       << m_nbUpdatableStates
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << m_nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  /// @todo only first order for now
  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER1 << "\n";

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(m_elementType[k].getNbNodesPerCell(),
				      CFPolyOrder::ORDER1,
				      m_dimension) << "\n";

  }

  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbCellsPerType() << " ";
  }
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    CFuint nbCellsPerType = m_elementType[k].getNbCellsPerType();
    CFuint nbNodesPerCell = m_elementType[k].getNbNodesPerCell();
    CFuint nbStatesPerCell  = nbNodesPerCell;

    for (CFuint i = 0; i < nbCellsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerCell; ++j) {
        fout << m_elementType[k].getTableConnectivity()(i,j) <<" " ;
      }
      for (CFuint j = 0; j < nbStatesPerCell; ++j) {
        fout << m_elementType[k].getTableConnectivity()(i,j) <<" " ;
      }
      fout << "\n";
    }
  }
  
  CFLog(VERBOSE, "THOR2CFmeshConverter::writeContinuousElements() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << m_isWithSolution << "\n";

  if (m_isWithSolution) {
    fout.precision(14);
    fout << *m_variables;
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_NODE " << "\n";
  for (CFuint k = 0; k < m_nbUpdatableNodes; ++k) {
    for (CFuint j = 0; j < m_dimension; ++j)
    {
      fout.precision(14);
      fout.setf(ios::scientific,ios::floatfield);
      fout << (*m_coordinate)(k,j) << " ";
    }
    fout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "THOR2CFmeshConverter::writeDiscontinuousElements() => start\n");
  
  // this is here for backward compatibility
  // these numbers are always zero
  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  // count the number of states and cells
  std::valarray<CFuint> columnPattern(m_nbCells);
  CFuint totalNbStates = 0;
  CFuint totalNbCells = 0;
  CFuint cellitr = 0;
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    const CFuint nbCellsPerType = m_elementType[k].getNbCellsPerType();
    const CFuint nbStatesPerCell = getNbStatesInType(k);

    totalNbCells += nbCellsPerType;
    totalNbStates += nbCellsPerType*nbStatesPerCell;
    for (; cellitr < totalNbCells; ++cellitr) columnPattern[cellitr] = nbStatesPerCell;
  }

  cf_assert(m_nbCells == totalNbCells);

  fout << "!NB_NODES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  fout << "!NB_STATES "
       << totalNbStates
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << m_nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) m_solOrder << "\n";

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(m_elementType[k].getNbNodesPerCell(),
				      CFPolyOrder::ORDER1,
				      m_dimension) << "\n";
  }

  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbCellsPerType() << " ";
  }
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    fout << m_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    fout << getNbStatesInType(k) << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";

  // create the connectivity cell to pair of node and state
  // for the discontinuous CFPolyOrder::ORDER1 case because the
  // the id of nodes is not the same as id of states
  if (m_solOrder != CFPolyOrder::ORDER0)
  {
    m_scon = new Common::Table< std::pair<CFuint,CFuint> > (columnPattern);
  }

  CFuint countElem = 0;
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    const CFuint nbCellsPerType = m_elementType[k].getNbCellsPerType();
    const CFuint nbNodesPerCell = m_elementType[k].getNbNodesPerCell();
    const CFuint nbStatesPerCell = getNbStatesInType(k);

    for (CFuint i = 0; i < nbCellsPerType; ++i)
    {
      for (CFuint j = 0; j < nbNodesPerCell; ++j)
      {
        fout << m_elementType[k].getTableConnectivity()(i,j) << " " ;
      }

      for (CFuint j = 0; j < nbStatesPerCell; ++j)
      {
        fout << countElem << " " ;
        // fill the connectivity because we need to write the faces
        if (m_solOrder != CFPolyOrder::ORDER0)
        {
          cf_assert(m_scon != CFNULL);
          (*m_scon)(i,j) = std::make_pair(m_elementType[k].getTableConnectivity()(i,j),countElem);
        }
        ++countElem;
      }
      fout << "\n";

    } // loop cells
  } // loop element types
  
  CFLog(VERBOSE, "THOR2CFmeshConverter::writeDiscontinuousElements() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

CFuint THOR2CFmeshConverter::getNbStatesInType(const CFuint& typeID)
{
 if(isDiscontinuous()) {
  if (m_solOrder == CFPolyOrder::ORDER0)
  {
    return 1;
  }
  else // m_solOrder == CFPolyOrder::ORDER1
  {
    return m_elementType[typeID].getNbNodesPerCell();
  }
 }
 else return m_elementType[typeID].getNbNodesPerCell();
}

//////////////////////////////////////////////////////////////////////////////

void THOR2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << m_isWithSolution << "\n";

  CFuint nbVariables = m_nbLamVariables + m_nbTurbVariables;
  if (nbVariables == 0)
  {
    nbVariables = PhysicalModelStack::getActive()->getNbEq();
  }
  std::valarray<CFreal> averageState(nbVariables);

  if (m_isWithSolution)
  {
    if (m_solOrder == CFPolyOrder::ORDER0) // one state per cell
    {
      for (CFuint k = 0; k < getNbElementTypes(); ++k)
      {
        CFuint nbCellsPerType = m_elementType[k].getNbCellsPerType();
        CFuint nbNodesPerCell = m_elementType[k].getNbNodesPerCell();

        cf_assert(nbCellsPerType > 0);
        for (CFuint i = 0; i < nbCellsPerType; ++i) {
          averageState = 0.0;
          for (CFuint j = 0; j < nbNodesPerCell; ++j) {
            const CFuint nodeID = m_elementType[k].getTableConnectivity()(i,j);
            const vector<CFreal>& nodalState = m_variables->getRow(nodeID);
            for(CFuint v = 0; v < nbVariables; ++v) {
              averageState[v] += nodalState[v];
            }
          }
          averageState /= nbNodesPerCell; // compute the average state
          for (CFuint iVar = 0; iVar < nbVariables; ++iVar) {
            fout.precision(14);
            fout.setf(ios::scientific,ios::floatfield);
            fout << averageState[iVar] << " ";
          }
          fout << "\n";
        }
      }
    }
    else // each state is on each node
    {
      for (CFuint k = 0; k < getNbElementTypes(); ++k)
      {
        CFuint nbCellsPerType = m_elementType[k].getNbCellsPerType();
        CFuint nbNodesPerCell = m_elementType[k].getNbNodesPerCell();
        cf_assert(nbNodesPerCell == getNbStatesInType(k));
        cf_assert(nbCellsPerType > 0);

        for (CFuint i = 0; i < nbCellsPerType; ++i) {
          for (CFuint j = 0; j < nbNodesPerCell; ++j) {
            const CFuint nodeID = m_elementType[k].getTableConnectivity()(i,j);
            const vector<CFreal>& nodalState = m_variables->getRow(nodeID);

            for (CFuint iVar = 0; iVar < nbVariables; ++iVar)
            {
              fout.precision(14);
              fout << nodalState[iVar] << " ";
            }
            fout << "\n";
        }
      }
    }
  }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace THOR2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
