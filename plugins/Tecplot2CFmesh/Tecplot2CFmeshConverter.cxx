// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>

#include "Common/CFLog.hh"
#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Common/Stopwatch.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/CFMultiMap.hh"
#include "Common/OSystem.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataStorage.hh"
#include "Tecplot2CFmesh/Tecplot2CFmeshConverter.hh"
#include "Tecplot2CFmesh/Tecplot2CFmesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Tecplot2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Tecplot2CFmeshConverter,
			    MeshFormatConverter,
			    Tecplot2CFmeshModule,
			    1>
tecplot2CFmeshConverterProvider("Tecplot2CFmesh");

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::defineConfigOptions(Config::OptionList& options)
{  
  options.addConfigOption< bool >
    ("HasBlockFormat","Data are read in blocks (1 block <=> 1 variable).");
  
  options.addConfigOption< vector<string> >
    ("SurfaceTRS","Names of the surface TRS.");
  
  options.addConfigOption< vector<string> >
    ("ReadVariables","List of the names of the variables to read.");
  
  options.addConfigOption< CFuint >
    ("NbElementTypes","Number of element types.");
  
  options.addConfigOption< CFuint >
    ("Precision","Number of digits to be considered for nodal coordinates matching.");
  
  options.addConfigOption< bool >
    ("HasAllSurfFile","Flag telling if a TECPLOT file *.allsurf.plt including all boundaries is given.");
  
  options.addConfigOption< bool >
    ("SaveNodalStates","Store the nodal states while converting the mesh.");
  
  options.addConfigOption< std::string >("InterpolateFrom","The donor file for interpolation");
}
      
//////////////////////////////////////////////////////////////////////////////
      
Tecplot2CFmeshConverter::Tecplot2CFmeshConverter (const std::string& name)
  : MeshFormatConverter(name),
    m_dimension(0),
    m_elemMapCoordToNodeID(),
    m_nodalVarPtr(),
    m_elementType(),
    m_hasBlockFormat(false)
{
  addConfigOptionsTo(this);
  m_hasBlockFormat = false;
  setParameter("HasBlockFormat", &m_hasBlockFormat);
  
  m_surfaceTRS = vector<string>();
  setParameter("SurfaceTRS", &m_surfaceTRS);
  
  m_readVars = vector<string>();
  setParameter("ReadVariables", &m_readVars);
  
  m_nbElemTypes = 1;
  setParameter("NbElementTypes", &m_nbElemTypes);
  
  m_precision = 8;
  setParameter("Precision", &m_precision); 
  
  m_hasAllSurfFile = false;
  setParameter("HasAllSurfFile", &m_hasAllSurfFile); 
  
  m_saveNodalStates = false;
  setParameter("SaveNodalStates", &m_saveNodalStates); 
  
  m_interpolateFromFileName = "";
  setParameter("InterpolateFrom", &m_interpolateFromFileName);
}
      
//////////////////////////////////////////////////////////////////////////////

Tecplot2CFmeshConverter::~Tecplot2CFmeshConverter()
{
  for (CFuint i = 0; i < m_elementType.size(); ++i) {
    deletePtr(m_elementType[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::configure ( Config::ConfigArgs& args )
{
  MeshFormatConverter::configure(args);  
  
  if (m_interpolateFromFileName.empty() && !m_hasAllSurfFile)
  {
    CFLog(WARN, "Tecplot2CFmeshConverter::configure() => no interpolation is requested\n");
  }
  
  // enforce flag to be true if a donor file for interpolation is specified
  if (!m_interpolateFromFileName.empty()) {m_hasAllSurfFile = true;}
  
  boost::filesystem::path fp (m_interpolateFromFileName);
  m_interpolateFromFileName = fp.string();
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::checkFormat(const boost::filesystem::path& filepath)
{
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::readTecplotFile(CFuint nbZones, 
					      const boost::filesystem::path& filepath, 
					      const string& extension, 
					      bool isBoundary) 
{
  using namespace boost::filesystem;
    
  // read the cells
  string line = "";
  CFuint countl = 0;
  vector<string> words;
  SelfRegistPtr<FileHandlerInput> fhandle = SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
  path meshFile = change_extension(filepath, extension);
  
  ifstream& fin = fhandle->open(meshFile);
  getTecplotWordsFromLine(fin, line, countl, words); // TITLE
  
  vector<string> vars;
  readVariables(fin, line, countl, words, vars); // VARIABLES
  for (CFuint i = 0; i < nbZones; ++i) {
    CFLog(VERBOSE, "START Reading ZONE " << i << "\n");
    CFLog(VERBOSE, "isBoundary =" << isBoundary << "\n");
    readZone(fin, countl, vars, isBoundary);
    CFLog(VERBOSE, "END   Reading ZONE " << i << "\n");
  }
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::readZone(ifstream& fin, 
				       CFuint& countl, 
				       vector<string>& vars,
				       bool isBoundary)
{  
  string line = "";
  vector<string> words;
  CFuint nbElems = 0; 
  CFuint nbNodes = 0; 
  string cellType = "";
  
  // keep on reading lines till when you get a line != empty
  while (line.find("ZONE") == string::npos) {
    CFLog(VERBOSE, "line = " << line << "\n"); 
    getTecplotWordsFromLine(fin, line, countl, words); // ZONE Description
  }
  CFLog(VERBOSE, CFPrintContainer<vector<string> >("line = ", &words) << "\n");

  bool foundN = false;
  bool foundE = false;
  bool foundET = false;
  
  while (!foundN || !foundE || !foundET) {
    string tmpStr = "";
    for (CFuint i = 0; i < words.size(); ++i) {
      CFLog(VERBOSE, "Tecplot2CFmeshConverter::readZone() => current word is [" << words[i] << "]\n");
      
      string nextWord = (i < words.size()-1) ? words[i+1] : "";
      if (words[i].find("N=") != string::npos) {
	getValueString(string("N="), words[i], nextWord, tmpStr);
	const CFuint nbN = StringOps::from_str<CFuint>(tmpStr);
	nbNodes = (nbN > 0) ? nbN : nbNodes;
	foundN = true;
      }
      else if (words[i].find("Nodes=") != string::npos) {
	getValueString(string("Nodes="), words[i], nextWord, tmpStr);
	const CFuint nbN = StringOps::from_str<CFuint>(tmpStr);
	nbNodes = (nbN > 0) ? nbN : nbNodes;
	foundN = true;
      }
      else if (words[i].find("NODES=") != string::npos) {
	getValueString(string("NODES="), words[i], nextWord, tmpStr);
        const CFuint nbN = StringOps::from_str<CFuint>(tmpStr);
        nbNodes = (nbN > 0) ? nbN : nbNodes;
        foundN = true;
      }
      else if (words[i].find("ZONETYPE=") != string::npos) {
	getValueString(string("ZONETYPE="), words[i], nextWord, cellType);
	foundET = true;
      }
      else if (words[i].find("E=") != string::npos) {
	getValueString(string("E="), words[i], nextWord, tmpStr);
	const CFuint nbC = StringOps::from_str<CFuint>(tmpStr);
	nbElems = (nbC > 0) ? nbC : nbElems;
	foundE = true;
      }
      else if (words[i].find("Elements=") != string::npos) {
	getValueString(string("Elements="), words[i], nextWord, tmpStr);
	const CFuint nbC = StringOps::from_str<CFuint>(tmpStr);
	nbElems = (nbC > 0) ? nbC : nbElems;
	foundE = true;
      }
      else if (words[i].find("ELEMENTS=") != string::npos) {
        getValueString(string("ELEMENTS="), words[i], nextWord, tmpStr);
        const CFuint nbC = StringOps::from_str<CFuint>(tmpStr);
        nbElems = (nbC > 0) ? nbC : nbElems;
        foundE = true;
      }
      else if (words[i].find("ET=") != string::npos) {
	getValueString(string("ET="), words[i], nextWord, cellType);
	foundET = true;
      }
      
      if (foundN && foundE && foundET) break; 
    } 
    
    if (foundN && foundE && foundET) break; 
    getTecplotWordsFromLine(fin, line, countl, words); // ZONE Description
  }
  
  // at this point in the file, you could have
  // DATAPACKING=POINT
  // DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )
  // in that case, read it and jump to the useful data
  long startData = fin.tellg();
  getTecplotWordsFromLine(fin, line, countl, words);
  if (line.find("DATAPACKING=") != string::npos) {
    startData = fin.tellg();
    getTecplotWordsFromLine(fin, line, countl, words);
  }
  //  if the line contains "DT=" read it, otherwise go back because this is data to be stored 
  if (line.find("DT=") == string::npos) {
    fin.seekg(startData);
  }

  m_dimension = max(m_dimension, getDim(cellType));
  const CFuint nbNodesInCell = getNbNodesInCell(cellType); 
    
  CFLog(VERBOSE, "ZONE has: nbNodes = " << nbNodes << ", nbElems = " << nbElems 
	<< ", cellType = " << cellType << ", nbNodesInCell = " << nbNodesInCell << "\n");
  
  const CFuint nbVarsToStore = (skipSolution()) ? m_dimension : (m_dimension + m_readVars.size());
  
  CFLog(VERBOSE, "Vars to read = ");
  for (CFuint v = 0; v < m_readVars.size(); ++v) {
    CFLog(VERBOSE, "(" << v << ")" << m_readVars[v] <<  " ");
  }
  CFLog(VERBOSE, "\n");
  
  ElementTypeTecplot* elements = new ElementTypeTecplot(nbElems, nbNodesInCell, nbNodes, nbVarsToStore);
  
  CFLog(VERBOSE, "m_dimension = " << m_dimension << ", m_readVars.size() = " 
	<< m_readVars.size() << ", nbVarsToStore = " << nbVarsToStore << "\n");	  
  
  // mesh/solution variables
  const CFuint nbVars = vars.size(); 
  CFreal value = 0;
  
  vector<CFint> varID(nbVars, -1);
  // AL: the first m_dimension variables are assumed to be ALWAYS coordinates x,y,(z)
  for (CFuint i = 0; i < m_dimension; ++i) {
    varID[i] = i;
  }
  
  NodeDim currNode(CFNULL);
  NodeDim::setSize(m_dimension);
  
  // find matching variable names and assign the corresponding variable IDs
  // but not for boundary zones
  if (!isBoundary) {
    for (CFuint i = m_dimension; i < nbVars; ++i) {
      for (CFuint iVar = 0; iVar < m_readVars.size(); ++iVar) {
	if (vars[i] == m_readVars[iVar]) {
	  varID[i] = iVar + m_dimension;
	}
      }
    }
  }
  
  cf_assert(nbNodes > 0);
  if (!isBoundary) {
    const CFuint newSize = m_nodalVarPtr.size() + nbNodes;
    m_nodalVarPtr.reserve(newSize);
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::readZone() => m_nodalVarPtr.capacity() = " << m_nodalVarPtr.capacity()  << "\n");
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::readZone() => m_nodalVarPtr.size() = " << m_nodalVarPtr.size()  << "\n");
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::readZone() => nbNodes = " << nbNodes << "\n");
    cf_assert(m_nodalVarPtr.capacity() == newSize);
    cf_assert(m_nodalVarPtr.capacity() > 0);
  }
  
  vector<CFuint> mapNodeID2UniqueNodeID(nbNodes);
  // const CFreal eps = std::pow(10.,-1.*m_precision);
  
  if (m_hasBlockFormat) {
    // first read all variables in each node
    // starting with all x's, all y's ...
    for (CFuint i = 0; i < nbVars; ++i) {
      for (CFuint n = 0; n < nbNodes; ++n) {
	fin >> value;
	if (varID[i] > -1) {
	  elements->setVar(n,varID[i], value);
	}
      }
    }
    
    // then map individual nodal coordinates (x,y,(z)) to node IDs
    if (!isBoundary) { 
      for (CFuint i = 0; i < nbNodes; ++i) {
	// add new entry into m_nodalVarPtr if the current nodal x,y,(z) are new
	CFreal *const nodalPtr = elements->getNodalVarPtr(i);
	currNode.reset(nodalPtr);
	if (m_elemMapCoordToNodeID.count(currNode) == 0) {
	  const CFuint maxNodeID = m_nodalVarPtr.size();
	  m_nodalVarPtr.push_back(nodalPtr);
	  
	  // reduce precision to allow for sufficiently high matching tolerance
	  // CFLog(INFO, "Rounding up node = " << currNode << " -> ");
	  // MathTools::MathFunctions::roundUpPrecision(currNode, eps);
	  // CFLog(INFO, currNode  << "\n");
	  m_elemMapCoordToNodeID.insert(pair<NodeDim,CFuint>(currNode, maxNodeID));
	  
	  mapNodeID2UniqueNodeID[i] = maxNodeID+1;
	  cf_assert(m_nodalVarPtr.size() == maxNodeID+1);
	}
	else {
	  cf_assert(m_elemMapCoordToNodeID.count(currNode) == 1);
	  mapNodeID2UniqueNodeID[i] = m_elemMapCoordToNodeID.find(currNode)->second + 1;
	}
      }
    }
  }
  else {
    for (CFuint i = 0; i < nbNodes; ++i) {
      if (!isBoundary) { 
	CFLog(DEBUG_MAX, "Node["<< i << "] => ");
	for (CFuint n = 0; n < nbVars; ++n) {
	  fin >> value;
	  CFLog(DEBUG_MAX, "(" << n << ": ");
	  if (varID[n] > -1) {
	    CFLog(DEBUG_MAX, value);
	    elements->setVar(i,varID[n], value);
	  }
	  CFLog(DEBUG_MAX, ") ");
	}
	CFLog(DEBUG_MAX, "\n");
	
	// add new entry into m_nodalVarPtr if the current nodal x,y,(z) are new
	CFreal *const nodalPtr = elements->getNodalVarPtr(i);
	currNode.reset(nodalPtr);
	CFLog(DEBUG_MAX, "Node[" << i << "] = " << currNode << " => ");
	// here we make sure not to insert twice the same node (as coming 
	// from different element zones sharing that node)
	if (m_elemMapCoordToNodeID.count(currNode) == 0) {
	  const CFuint maxNodeID = m_nodalVarPtr.size();
	  m_nodalVarPtr.push_back(nodalPtr);
	  m_elemMapCoordToNodeID.insert(pair<NodeDim,CFuint>(currNode, maxNodeID));
	  mapNodeID2UniqueNodeID[i] = maxNodeID+1;
	  CFLog(DEBUG_MAX, " IN  " << mapNodeID2UniqueNodeID[i] << "\n\n");
 	  cf_assert(m_nodalVarPtr.size() == maxNodeID+1);
	}
	else {
	  cf_assert(m_elemMapCoordToNodeID.count(currNode) == 1);
	  mapNodeID2UniqueNodeID[i] = m_elemMapCoordToNodeID.find(currNode)->second + 1;
	  CFLog(DEBUG_MAX, " OUT " << mapNodeID2UniqueNodeID[i] << "\n\n");
	}
      }
      else {
	for (CFuint n = 0; n < nbVars; ++n) { 
	  fin >> value;
	  if (n < m_dimension && varID[n] > -1) { 
	    elements->setVar(i,varID[n], value); 
	  } 
	}
      }
    }
  }
  
  // IMPORTANT: at this point TECPLOT numbering is used (+1 compared to CF numbering)
  // elements connectivity 
  if (getNbElementTypes() == 1) {
    for (CFuint i = 0; i < nbElems; ++i) {
      for (CFuint n = 0; n < nbNodesInCell; ++n) {
	fin >> (*elements)(i,n);
	elements->setBkpElementNode(i, n, (*elements)(i,n));
	if ((*elements)(i,n) > elements->getNbNodes()) {
	  CFLog(WARN, (*elements)(i,n) << " > " << elements->getNbNodes() << "\n");
	  cf_assert((*elements)(i,n) < elements->getNbNodes());
	}
      }
    }
  }
  else {
    for (CFuint i = 0; i < nbElems; ++i) {
      for (CFuint n = 0; n < nbNodesInCell; ++n) {
	fin >> (*elements)(i,n);
	elements->setBkpElementNode(i, n, (*elements)(i,n));
	cf_assert((*elements)(i,n) <= nbNodes); // TECPLOT node numbering starts from 1
	
	if (!isBoundary) { 
	  // reset the element node ID, considering only unique nodes 
	  elements->setNodeDim(n, currNode);
	  cf_assert(m_elemMapCoordToNodeID.count(currNode) == 1);
	  // node ID for Tecplot start from +1 (not 0)
	  (*elements)(i,n) = mapNodeID2UniqueNodeID[(*elements)(i,n)-1];
	  
	  // here the numbering is still TECPLOT's
	  if ((*elements)(i,n) > m_nodalVarPtr.size()) {
	    CFLog(WARN, (*elements)(i,n) << " > " << m_nodalVarPtr.size() << "\n");
	    cf_assert((*elements)(i,n) < m_nodalVarPtr.size());
	  }
	}
	else {
	  // on the boundary, we only check against max number of nodes in the zone 
	  if ((*elements)(i,n) > elements->getNbNodes()) {
	    CFLog(WARN, (*elements)(i,n) << " > " << elements->getNbNodes() << "\n");
	    cf_assert((*elements)(i,n) < elements->getNbNodes());
	  }
	}
      }
    }
  }
  
  // update the number of non-duplicated nodes per element
  elements->setNbUniqueNodesPerElem();
  
  m_elementType.push_back(elements);
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::convertBack(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  
  //  only reads if not yet been read
  //  readFiles(filepath);
  //  writeTecplot(filepath);
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeTecplot(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

//   using namespace boost::filesystem;
//   path outFile = change_extension(filepath, getOriginExtension());
  
//   Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
//   ofstream& fout = fhandle->open(outFile);
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeContinuousTrsData() => START\n");
  
  renumberTRSData();
  
  const CFuint startTRS = (m_hasAllSurfFile) ? m_nbElemTypes+1 : m_nbElemTypes;
  const CFuint nbTRSs = m_elementType.size() - startTRS;
  fout << "!NB_TRSs " << nbTRSs << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    const CFuint nbTRsInTRS = 1;
    fout << "!TRS_NAME " << m_surfaceTRS[iTRS] << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";
    
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
    const CFuint nbTRSFaces = trs->getNbElems();
    fout << nbTRSFaces << " ";
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";
    const CFuint nbNodesPerFace = trs->getNbNodesPerElem();
    for (CFuint iFace = 0; iFace < nbTRSFaces; ++iFace) {
      fout << nbNodesPerFace << " " << nbNodesPerFace << " ";
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	fout << (*trs)(iFace,iNode) << " ";
      }
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	fout << (*trs)(iFace,iNode) << " ";
      }
      fout << "\n";
    }
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeContinuousTrsData() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
{  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousTrsData() => START\n");
  
  renumberTRSData();
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousTrsData() 1\n");
  
  buildBFaceNeighbors();
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousTrsData() 2\n");
  
  const CFuint startTRS = (m_hasAllSurfFile) ? m_nbElemTypes+1 : m_nbElemTypes;
  const CFuint nbTRSs = m_elementType.size() - startTRS;
  fout << "!NB_TRSs " << nbTRSs << "\n";
  
  CFuint maxNbNodesInFace = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
    maxNbNodesInFace = std::max(maxNbNodesInFace, trs->getNbNodesPerElem());
  }
  
  set<CFuint> nodeIDList;
  vector<CFuint> localElemNodeIdx;
  localElemNodeIdx.reserve(maxNbNodesInFace);
  
  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    const CFuint nbTRsInTRS = 1;
    fout << "!TRS_NAME " << m_surfaceTRS[iTRS] << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";
    
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
    const CFuint nbTRSFaces = trs->getNbElems();
    fout << nbTRSFaces << " ";
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";
    
    // the same TRS can have different number of nodes per face
    for (CFuint iFace = 0; iFace < nbTRSFaces; ++iFace) {
      trs->getUniqueNodesInSingleElem(iFace, nodeIDList, localElemNodeIdx);
      const CFuint nbNodesPerFace = localElemNodeIdx.size();
      cf_assert(nbNodesPerFace > 0);
      cf_assert(nbNodesPerFace <= trs->getNbNodesPerElem());
      fout << nbNodesPerFace << " " << 1 << " ";
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	const CFuint localNodeID = localElemNodeIdx[iNode];
	fout << (*trs)(iFace,localNodeID) << " ";	
      }
      fout << trs->getNeighborID(iFace) << "\n";
    }
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousTrsData() => END\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::adjustToCFmeshNodeNumbering() 
{
  offsetNumbering(-1);
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::offsetNumbering(CFint offset)
{
  for (CFuint i = 0; i < m_elementType.size(); ++i) {
    m_elementType[i]->offsetNumbering(offset);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::readFiles() => START\n");
  
  try {
    CFLog(VERBOSE, CFPrintContainer<vector<string> >("surfaces to be read = ", &m_surfaceTRS) <<  "\n");
    
    // if a donor .plt file is specified, then a Tecplot macro interpolates the solution
    // and writes out all the boundaries of the interpolated mesh in .allsurf.plt file 
    if (!m_interpolateFromFileName.empty()) {
      interpolateTecplotSolution(filepath);
    }
    
    readTecplotFile(m_nbElemTypes, filepath, getOriginExtension(), false);
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::readFiles() => nb element types = " << m_elementType.size() << "\n");
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::readFiles() => nb unique nodes  = " << m_nodalVarPtr.size() << "\n");
    
    if (m_hasAllSurfFile) {
      CFLog(VERBOSE, "Tecplot2CFmeshConverter::readFiles() => reading .allsurf.plt\n");
      readTecplotFile(1, filepath, ".allsurf.plt", true);
    }
    
    readTecplotFile(m_surfaceTRS.size(), filepath, ".surf.plt", true);
  }
  catch (Common::Exception& e) {
    CFout << e.what() << "\n";
    throw;
  }
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::readFiles() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeContinuousElements(ofstream& fout)
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeContinuousElements() => START\n");
  
  const CFuint nbElemTypes = getNbElementTypes();
  CFuint nbElems = 0;
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    nbElems += m_elementType[k]->getNbElems();  
  } 
  cf_assert(nbElems > 0);
  
  const CFuint nbNodes = m_nodalVarPtr.size();
  cf_assert(nbNodes > 0);
  if (nbElemTypes == 1) {
    cf_assert(nbNodes == m_elementType[0]->getNbNodes());
  }
  
  fout << "!NB_NODES " << nbNodes << " " << 0 << "\n";
  fout << "!NB_STATES " << nbNodes << " " << 0 << "\n";
  fout << "!NB_ELEM "   << nbElems << "\n";
  fout << "!NB_ELEM_TYPES " << getNbElementTypes() << "\n";
  
  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER1 << "\n";
  
  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(m_elementType[k]->getNbNodesPerElem(),
				      CFPolyOrder::ORDER1,
				      m_dimension) << "\n";
  }
  
  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k]->getNbElems() << " ";
  }
  fout << "\n";
  
  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k]->getNbNodesPerElem() << " ";
  }
  fout << "\n";
  
  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k]->getNbNodesPerElem() << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    const CFuint nbElemsPerType  = m_elementType[k]->getNbElems();
    const CFuint nbNodesPerElem  = m_elementType[k]->getNbNodesPerElem();
    const CFuint nbStatesPerElem = nbNodesPerElem;
    
    // decrease IDs by 1 because Tecplot numbering starts from 1  
    for (CFuint i = 0; i < nbElemsPerType; ++i) {
      for (CFuint j = 0; j < nbNodesPerElem; ++j) {
        fout << (*m_elementType[k])(i,j) <<" " ;
      }
      for (CFuint j = 0; j < nbStatesPerElem; ++j) {
        fout << (*m_elementType[k])(i,j) <<" " ;
      }
      fout << "\n";
    }
  } 
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeContinuousElements() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeContinuousStates() => START\n");
  
  fout << "!LIST_STATE " << !skipSolution() << "\n";
   
  if (!skipSolution()) {
    const CFuint nbNodes = m_nodalVarPtr.size();
    const CFuint nbVars = m_elementType[0]->getNbVars();
    const CFuint startVar = m_dimension;
    for (CFuint k = 0; k < nbNodes; ++k) {
      cf_assert(k < m_nodalVarPtr.size());
      const CFreal *const startCoord = m_nodalVarPtr[k];
      for (CFuint j = startVar; j < nbVars; ++j) {
       	fout.precision(14);
	fout.setf(ios::scientific,ios::floatfield);
	fout << startCoord[j] << " ";
      }
      fout << "\n";
    }
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeContinuousStates() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void Tecplot2CFmeshConverter::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeNodes() => START\n");
  
  fout << "!LIST_NODE " << "\n";
  
  const CFuint nbNodes = m_nodalVarPtr.size();
  for (CFuint k = 0; k < nbNodes; ++k) {
    cf_assert(k < m_nodalVarPtr.size());
    const CFreal *const startCoord = m_nodalVarPtr[k];
    for (CFuint j = 0; j < m_dimension; ++j) {
      fout.precision(14);
     fout.setf(ios::scientific,ios::floatfield);
     fout << startCoord[j] << " ";
    }
    fout << "\n";
  }
  
  const CFuint nbVars = m_elementType[0]->getNbVars();
  const CFuint startVar = m_dimension;
  const CFuint nbNodalVars = nbVars - startVar; 
  
  // AL: here we store the nodal variables assuming that the global node IDs
  // won't change later (e.g. sue to mesh partitioning or renumbering) otherwise 
  // the DataHandle "initialNodalStates" will also need some ID mapping to be accessible
  if (m_saveNodalStates && (!skipSolution()) && nbNodalVars > 0) {
    const string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
    const string dhName = nsp + "_initialNodalStates";
    DataHandle<CFreal> initialNodalStates =
      MeshDataStack::getActive()->getDataStorage()->createData<CFreal>
      (dhName, nbNodes*nbNodalVars); 
    
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeNodes() => " << 
	  "[startVar, nbVars, nbNodalVars] = [" << startVar << ", " << nbVars << ", " << nbNodalVars << "]\n");
    
    for (CFuint k = 0; k < nbNodes; ++k) {
      cf_assert(k < m_nodalVarPtr.size());
      const CFreal *const startCoord = m_nodalVarPtr[k];
      const CFuint firstVar = k*nbNodalVars;
      for (CFuint j = startVar; j < nbVars; ++j) {
	const CFuint idx = firstVar + j - startVar;
	cf_assert(idx < initialNodalStates.size());
	initialNodalStates[idx] = startCoord[j];
      }
    }
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeNodes() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousElements() => START\n");
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousElements() => total nb zones = " <<
	m_elementType.size() << "\n");
  // nb zones includes cells and boundary faces zones
  cf_assert(m_elementType.size() > m_nbElemTypes);
  
  const CFuint nbElemTypes = getNbElementTypes();
  CFuint nbElems = 0;
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    nbElems += m_elementType[k]->getNbElems();  
  }
  cf_assert(nbElems > 0);
  
  const CFuint nbNodes = m_nodalVarPtr.size();
  cf_assert(nbNodes > 0);
  if (nbElemTypes == 1) {
    cf_assert(nbNodes == m_elementType[0]->getNbNodes());
  }
  
  fout << "!NB_NODES " << nbNodes << " " << 0 << "\n";
  fout << "!NB_STATES " << nbElems << " " << 0 << "\n";
  fout << "!NB_ELEM "   << nbElems << "\n";
  fout << "!NB_ELEM_TYPES " << nbElemTypes << "\n";
  
  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER0 << "\n";
  
  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    fout << MapGeoEnt::identifyGeoEnt(m_elementType[k]->getNbUniqueNodesPerElem(),
				      CFPolyOrder::ORDER1,
				      m_dimension) << " ";
  }
  fout << "\n";
  
  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    fout << m_elementType[k]->getNbElems() << " ";
  }
  fout << "\n";
  
  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    fout << m_elementType[k]->getNbUniqueNodesPerElem() << " ";
  }
  fout << "\n";
  
  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    fout << 1 << " ";
  }
  fout << "\n";
  
  fout << "!LIST_ELEM " << "\n";

  CFuint elemOffset = 0;
  for (CFuint k = 0; k < nbElemTypes; ++k) {
    const CFuint nbElemsPerType  = m_elementType[k]->getNbElems();
    const CFuint nbNodesPerElem  = m_elementType[k]->getNbNodesPerElem(); 
    const CFuint nbUniqueNodesPerElem  = m_elementType[k]->getNbUniqueNodesPerElem();
    const CFuint nbStatesPerElem = 1;
    // decrease IDs by 1 because Tecplot numbering starts from 1 
    if (nbNodesPerElem == nbUniqueNodesPerElem) {
      for (CFuint i = 0; i < nbElemsPerType; ++i) {
	for (CFuint j = 0; j < nbNodesPerElem; ++j) {
	  fout << (*m_elementType[k])(i,j) <<" " ;
	}
	for (CFuint j = 0; j < nbStatesPerElem; ++j) {
	  fout << i+elemOffset <<" " ;
	}
	fout << "\n";
      }
    }
    else {
      // case of pyramid or prism: here filter out duplicated node IDs 
      cf_assert(nbNodesPerElem > nbUniqueNodesPerElem);
      const vector<CFuint>& localElemNodeIdx = m_elementType[k]->getLocalElemNodeIdx();
      for (CFuint i = 0; i < nbElemsPerType; ++i) {
	for (CFuint j = 0; j < nbUniqueNodesPerElem; ++j) {
	  const CFuint localNodeID = localElemNodeIdx[j];
	  fout << (*m_elementType[k])(i,localNodeID) <<" " ;
	}
	for (CFuint j = 0; j < nbStatesPerElem; ++j) {
	  fout << i+elemOffset <<" " ;
	}
	fout << "\n";
      }
    }
    elemOffset += nbElemsPerType;
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousElements() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousStates() => START\n");
  
  fout << "!LIST_STATE " << !skipSolution() << "\n";
  
  if (!skipSolution()) {
    cf_assert(m_readVars.size() > 0);
    RealVector averageState(m_readVars.size());
    
    const CFuint nbElemTypes = getNbElementTypes();
    CFuint offset = 0;
    for (CFuint e = 0; e < nbElemTypes; ++e) {
      SafePtr<ElementTypeTecplot> elements = m_elementType[e]; 
      
      // consider the first element (all the other in this type are the same)
      // keep track of the local (i.e. inside the element itself) node IDs which have
      // to be written avoiding to consider duplicates
      const vector<CFuint>& localElemNodeIdx = elements->getLocalElemNodeIdx();
      const CFuint nbElems = elements->getNbElems();
      const CFuint nbVars = elements->getNbVars();
      const CFuint startVar = m_dimension;
      const CFuint nbUniqueNodesPerElem = elements->getNbUniqueNodesPerElem();
      const CFreal invNbNodeInElem = 1./static_cast<CFreal>(nbUniqueNodesPerElem); 
      for (CFuint k = 0; k < nbElems; ++k) {
	// first compute the average state out of the nodal values
	averageState = 0.;
	for (CFuint iNode = 0; iNode < nbUniqueNodesPerElem; ++iNode) {
	  const CFuint localNodeID = localElemNodeIdx[iNode];
	  const CFuint nodeID = elements->getBkpElementNode(k, localNodeID)-1;
	  CFLog(DEBUG_MAX, e << " => (" << k << ", " << localNodeID << ") => " 
		<< nodeID << " < " << elements->getNbNodes() << "\n");
	  cf_assert(nodeID < elements->getNbNodes());
	  for (CFuint j = startVar; j < nbVars; ++j) {
	    averageState[j-m_dimension] += elements->getVar(nodeID,j);
	  }
	}
	averageState *= invNbNodeInElem;
	
	fout.precision(14);
	fout.setf(ios::scientific,ios::floatfield);
	fout << averageState << endl;
      }
    }
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::writeDiscontinuousStates() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

CFuint Tecplot2CFmeshConverter::getNbNodesInCell(const std::string& etype) const
{
  if (etype == "LINESEG" || etype == "FELineseg" || 
      etype == "FELineSeg" || etype == "FELINESEG") {
    return 2;
  }
  else if (etype == "TRIANGLE" || etype == "FETriangle" || etype == "FETRIANGLE" ) {
    return 3;
  }
  else if (etype == "QUADRILATERAL" || etype == "TETRAHEDRON" || 
	   etype == "FEQuadrilateral" || etype == "FETetrahedron" || 
	   etype == "FEQUADRILATERAL" || etype == "FETETRAHEDRON") {
    return 4;
  }
  //  else if (etype == "BRICK" || etype == "FEBrick") {
  // AL: note prysm, pyramids, FEPOLYGON, FEPOLYHEDRON not supported
  return 8;
}

//////////////////////////////////////////////////////////////////////////////

CFuint Tecplot2CFmeshConverter::getDim(const std::string& etype) const
{
  if (etype == "LINESEG" || etype == "FELineseg" || 
      etype == "FELineSeg" || etype == "FELINESEG") {
    return DIM_1D;
  }
  else if (etype == "TRIANGLE" || etype == "QUADRILATERAL" || 
	   etype == "FETriangle" || etype == "FEQuadrilateral" ||
	   etype == "FETRIANGLE" || etype == "FEQUADRILATERAL") {
    return DIM_2D;
  }
  return DIM_3D;
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::getValueString(const std::string& key, 
					     const std::string& in,
					     const std::string& next,
					     std::string& out) const
{
  // if the given key does not coincide with the full word
  const string commaKey = "," + key;
  if (key != in && commaKey != in) {
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::getValueString() => key = [" <<  key << "] != [" << in << "]\n");
    const size_t pos = in.find(key);
    if (pos != std::string::npos) {
      const size_t posComma = in.find(",");
      size_t length = 0;
      if (posComma > pos) {
	length = min(in.size()-key.size(), posComma-key.size());
      }
      else {
	// this tackles badly formatted cases, e.g. [F=FEPOINT,ET=TRIANGLE] where [ET=] is the key
	length = in.size()-pos-key.size();
      }
      out = in.substr(pos+key.size(), length);
    }
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::getValueString() => out = [" <<  out << "]\n");
  }
  else {
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::getValueString() => key = [" <<  key << "] = [" << in << "]\n");
    // if the given key coincides with the full word
    // the next word gives the value
    cf_assert(next != "");
    // if there is a "," it can only be at the end of the string
    // this should tackle cases like " ,N=   3000," or ", N= 3000 ," or ", N= 3000,"  
    const size_t length = (next.find(",") == std::string::npos) ? next.size() : next.size()-1;
    out = next.substr(0, length);
    CFLog(VERBOSE, "Tecplot2CFmeshConverter::getValueString() => out = [" <<  out << "]\n");
  }
}
      
//////////////////////////////////////////////////////////////////////////////
 
void Tecplot2CFmeshConverter::readVariables(ifstream& fin,
					    std::string& line,
					    CFuint&  lineNb,
					    vector<string>& words, 
					    vector<string>& vars)
{  
  long startZone = fin.tellg();
  getTecplotWordsFromLine(fin, line, lineNb, words); // VARIABLES
  
  while (line.find("ZONE") == string::npos) {
    // skip line if the keyword "FILETYPE" is found
    if (line.find("FILETYPE") == string::npos) {
      for (CFuint i = 0; i < words.size(); ++i) {
	if (words[i].find('=')==string::npos && words[i].find("VARIABLES")==string::npos) {
	  CFLog(VERBOSE, "detected variable named <" << words[i] << "> \n");
	  // store all the variable names
	  vars.push_back(words[i]);
	}
      }
      
      startZone = fin.tellg(); 
    }
    getTecplotWordsFromLine(fin, line, lineNb, words); // VARIABLES
  }
  CFLog(VERBOSE, "\n");
  
  // we are assuming that the first "ZONE" comes right after the list of variables
  // set ifstream pointer to the position right before "ZONE"
  fin.seekg(startZone);
}
      
//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::buildBFaceNeighbors()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::buildBFaceNeighbors() =>START\n");
  
  const CFuint nbNodes = m_nodalVarPtr.size();
  vector<bool> isBNode(nbNodes, false);
  
  // zonesID=[0,startTRS) are cells
  // boundary TRS data start at zoneID=startTRS
  const CFuint startTRS = (m_hasAllSurfFile) ? m_nbElemTypes+1 : m_nbElemTypes;
  const CFuint nbTRSs = m_elementType.size() - startTRS;
  
  //count the number of boundary faces
  CFuint nbBFaces = 0;
  
  // create mapping between boundary node IDs and faces (trsID, faceID)
  CFMultiMap<CFuint, pair<CFuint, CFuint> > mapBNode2BFace;
  typedef CFMultiMap<CFuint, pair<CFuint, CFuint> >::MapIterator MapIterator;
  
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
    const CFuint nbTrsFaces = trs->getNbElems(); 
    const CFuint nbNodesPerFace = trs->getNbNodesPerElem();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	// this nodeID must be the global node ID (a renumbering has been performed)
	const CFuint nodeID = (*trs)(iFace,iNode);
	cf_assert(nodeID < nbNodes);
	
	// the map stores the TRS ID (starting from the last cell TRS ID)  
	mapBNode2BFace.insert(nodeID, pair<CFuint,CFuint>(startTRS + iTRS,iFace));
	isBNode[nodeID] = true;
      }
    }
    nbBFaces += nbTrsFaces;
  }
  mapBNode2BFace.sortKeys();
  
  CFuint countNeighborIDs = 0;
  CFuint elemOffset = 0;
  for (CFuint k = 0; k < m_nbElemTypes; ++k) {
    SafePtr<ElementTypeTecplot> elements = m_elementType[k]; 
    // loop over elements
    const CFuint nbElems = elements->getNbElems(); 
    const CFuint nbNodesInElem = elements->getNbNodesPerElem(); 
    
    for (CFuint e = 0; e < nbElems; ++e) {
      for (CFuint n = 0; n < nbNodesInElem; ++n) {
	const CFuint nodeID = (*elements)(e,n);
	if (isBNode[nodeID]) {
	  bool foundNode = false;
	  pair<MapIterator, MapIterator> faces = mapBNode2BFace.find(nodeID, foundNode);
	  cf_assert(foundNode);
	  
	  if (foundNode) {
	    // loop over all faces containing nodeID
	    bool faceFound = false;
	    for (MapIterator faceInMapItr = faces.first; 
		 faceInMapItr != faces.second && (faceFound == false);
		 ++faceInMapItr) {
	      
	      const CFuint trsID  = faceInMapItr->second.first;
	      const CFuint faceID = faceInMapItr->second.second;
	      cf_assert(trsID > 0);
	      SafePtr<ElementTypeTecplot> currTRS = m_elementType[trsID];
	      
	      // proceed only if the neighbor hasn't been set yet
	      if (currTRS->getNeighborID(faceID) == -1) {
		const CFuint nbENodes =  currTRS->getNbNodesPerElem(); 
		CFuint count = 0;
		
		// check if all nodes of the current face belong to the corrent element
		for (CFuint in = 0; in < nbENodes; ++in) {
		  for (CFuint en = 0; en < nbNodesInElem; ++en) {
		    if ((*currTRS)(faceID, in) == (*elements)(e,en)) {
		      count++;
		      break;
		    }
		  }
		}
		
		if (count == nbENodes) {
		  faceFound = true;
		  currTRS->setNeighborID(faceID, e+elemOffset);
		  countNeighborIDs++;
		}
	      }
	    }
	  }
	}  
      }
    }
    elemOffset += nbElems;
  }
  
  if (countNeighborIDs != nbBFaces) {
    CFLog(ERROR, "countNeighborIDs != nbBFaces => " << countNeighborIDs << " != " << nbBFaces << "\n"); 
    exit(1);
  } 
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::buildBFaceNeighbors() => END\n"); 
}

//////////////////////////////////////////////////////////////////////////////

void Tecplot2CFmeshConverter::interpolateTecplotSolution(const boost::filesystem::path& filepath)
{
  using namespace boost::filesystem;
  
  CFLog(INFO, "Tecplot2CFmeshConverter::interpolateTecplotSolution() => START\n");
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  path meshFile = change_extension(filepath, getOriginExtension());
  
  SelfRegistPtr<FileHandlerOutput> fhandle = SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();
  
  // string macroFileName = "interpolate.mcr";
  // path macroFile = DirPaths::getInstance().getWorkingDir()/macroFileName;
  // ofstream& fout = fhandle->open(macroFile);
  ofstream fout("interpolate.mcr");
  
  fout << "#!MC 1200\n";
  fout << "# Created by Tecplot 360 build 12.0.0.3454\n";
  fout << "$!VarSet |MFBD| = '" << DirPaths::getInstance().getWorkingDir().string() << "'\n";  
  fout << "$!READDATASET  '\"|MFBD|/" << m_interpolateFromFileName << "\" '\n";
  fout << "  READDATAOPTION = NEW\n";
  fout << "  RESETSTYLE = YES\n";
  fout << "  INCLUDETEXT = NO\n";
  fout << "  INCLUDEGEOM = NO\n";
  fout << "  INCLUDECUSTOMLABELS = NO\n";
  fout << "  VARLOADMODE = BYNAME\n";
  fout << "  ASSIGNSTRANDIDS = YES\n";
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  if (dim == DIM_3D) {
    fout << "  INITIALPLOTTYPE = CARTESIAN3D\n";
  }
  else if (dim == DIM_2D) {
    fout << "  INITIALPLOTTYPE = CARTESIAN2D\n";
  }
  fout << "  VARNAMELIST = '";
  for (CFuint i = 0 ; i < dim; ++i) {
    fout << "\"x" << i << "\" ";
  }
  for (CFuint i = 0 ; i < m_readVars.size(); ++i) {
    fout << m_readVars[i];
    if (i < m_readVars.size()-1) fout << " ";
  }
  fout << "'\n";

#ifdef CF_HAVE_BOOST_1_60
  fout << "$!READDATASET  '\"|MFBD|/" << meshFile.filename().string() << "\" '\n";
#else
#ifdef CF_HAVE_BOOST_1_59
  fout << "$!READDATASET  '\"|MFBD|/" << meshFile.filename().string() << "\" '\n";
#else
#ifdef CF_HAVE_BOOST_1_55
  fout << "$!READDATASET  '\"|MFBD|/" << meshFile.filename().string() << "\" '\n";
#else
#ifdef CF_HAVE_BOOST_1_54
  fout << "$!READDATASET  '\"|MFBD|/" << meshFile.filename().string() << "\" '\n";
#else
  CFLog(ERROR, "Tecplot2CFmeshConverter::interpolateTecplotSolution() => you need BOOST version >= 1.54 for this!\n");
  exit(1);
#endif
#endif
#endif
#endif
  
  fout << "  READDATAOPTION = APPEND\n";
  fout << "  RESETSTYLE = NO\n";
  fout << "  INCLUDETEXT = NO\n";
  fout << "  INCLUDEGEOM = NO\n";
  fout << "  INCLUDECUSTOMLABELS = NO\n";
  fout << "  VARLOADMODE = BYNAME\n";
  fout << "  ASSIGNSTRANDIDS = YES\n";
  fout << "  INITIALPLOTTYPE = CARTESIAN" << dim << "D\n"; // AL" check this for 2D
  fout << "  VARNAMELIST = '";
  for (CFuint i = 0 ; i < dim; ++i) {
    fout << "\"x" << i << "\" ";
  }
  for (CFuint i = 0 ; i < m_readVars.size(); ++i) {
    fout << m_readVars[i];
    if (i < m_readVars.size()-1) fout << " ";
  }
  fout << "'\n";
  fout << "$!INVERSEDISTINTERPOLATE\n"; 
  fout << "  SOURCEZONES =  [1]\n";
  fout << "  DESTINATIONZONE = 2\n";
  fout << "  VARLIST =  [" << dim+1 << "-" << dim+m_readVars.size() << "]\n";
  fout << "  INVDISTEXPONENT = 3.5\n";
  fout << "  INVDISTMINRADIUS = 0\n";
  fout << "  INTERPPTSELECTION = OCTANTNPOINTS\n";
  fout << "  INTERPNPOINTS = 8\n";
  fout << "$!CREATEFEBOUNDARY\n";
  fout << "  SOURCEZONE = 2\n";
  fout << "  REMOVEBLANKEDSURFACES = NO\n";

#ifdef CF_HAVE_BOOST_1_60
  fout << "$!WRITEDATASET  \"|MFBD|/" << meshFile.filename().string() << "\"\n";
#else
#ifdef CF_HAVE_BOOST_1_59
  fout << "$!WRITEDATASET  \"|MFBD|/" << meshFile.filename().string() << "\"\n";
#else
#ifdef CF_HAVE_BOOST_1_55
  fout << "$!WRITEDATASET  \"|MFBD|/" << meshFile.filename().string() << "\"\n";
#else
#ifdef CF_HAVE_BOOST_1_54
  fout << "$!WRITEDATASET  \"|MFBD|/" << meshFile.filename().string() << "\"\n";
#else
  CFLog(ERROR, "Tecplot2CFmeshConverter::interpolateTecplotSolution() => you need BOOST version >= 1.54 for this!\n");
  exit(1);
#endif
#endif
#endif
#endif

  fout << "  INCLUDETEXT = NO\n";
  fout << "  INCLUDEGEOM = NO\n";
  fout << "  INCLUDECUSTOMLABELS = NO\n";
  fout << "  ASSOCIATELAYOUTWITHDATAFILE = NO\n";
  fout << "  ZONELIST =  [2]\n";
  fout << "  BINARY = NO\n";
  fout << "  USEPOINTFORMAT = YES\n";
  fout << "  PRECISION = 9\n";
  fout << "  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT\n";

  path allSurfFile = change_extension(filepath, "allsurf.plt");
 
#ifdef CF_HAVE_BOOST_1_60
  fout << "$!WRITEDATASET  \"|MFBD|/" << allSurfFile.filename().string() << "\"\n";
#else 
#ifdef CF_HAVE_BOOST_1_59
  fout << "$!WRITEDATASET  \"|MFBD|/" << allSurfFile.filename().string() << "\"\n";
#else
#ifdef CF_HAVE_BOOST_1_55
  fout << "$!WRITEDATASET  \"|MFBD|/" << allSurfFile.filename().string() << "\"\n";
#else
#ifdef CF_HAVE_BOOST_1_54
  fout << "$!WRITEDATASET  \"|MFBD|/" << allSurfFile.filename().string() << "\"\n";
#else
  CFLog(ERROR, "Tecplot2CFmeshConverter::interpolateTecplotSolution() => you need BOOST version >= 1.54 for this!\n");
  exit(1);
#endif
#endif
#endif
#endif

  fout << "  INCLUDETEXT = NO\n";
  fout << "  INCLUDEGEOM = NO\n";
  fout << "  INCLUDECUSTOMLABELS = NO\n";
  fout << "  ASSOCIATELAYOUTWITHDATAFILE = NO\n";
  fout << "  ZONELIST =  [3]\n";
  fout << "  BINARY = NO\n";
  fout << "  USEPOINTFORMAT = YES\n";
  fout << "  PRECISION = 9\n";
  fout << "  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT\n";
  fout << "$!RemoveVar |MFBD|\n";
  fout << "$!QUIT\n";
  fout.close();
  
  // this suppose to have access to tec360 executable
  const string batchcommand = "tec360 -b -p interpolate.mcr";
  Common::OSystem::getInstance().executeCommand(batchcommand);
  
  CFLog(INFO,"Tecplot2CFmeshConverter::interpolateTecplotSolution() took " << stp.read() << "s\n");
}
      
//////////////////////////////////////////////////////////////////////////////
     
void Tecplot2CFmeshConverter::renumberTRSData()
{
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::renumberTRSData() => START/n");
  
  // assume 1 element type
  const CFuint totNbNodes = m_nodalVarPtr.size();
  const CFuint nbElemTypes = getNbElementTypes();
  
  // if one single TRS coming from TECPLOT interpolation is available, loop first 
  // on all boundary nodes to detect the global node IDs and be sure that the 
  // coordinates match between interior and boundary TRS
  // otherwise each TRS is processed, one at a time 
  vector<CFuint> allTRSNodeIDs;
  vector<CFreal> allTRSNodes; 
  CFuint counter = 0;
  CFuint totalNbTRSNodes = 0;
  const CFuint startTRS = m_nbElemTypes;
  const CFuint nbTRSs = (m_hasAllSurfFile) ? 1 : m_elementType.size() - startTRS;
  
  if (m_hasAllSurfFile) {
    cf_assert(nbElemTypes == 1);
    // preallocation of array to cache all node coordinates
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS];
    const CFuint nbTRSFaces = trs->getNbElems();
    const CFuint nbNodesPerFace = trs->getNbNodesPerElem();
    allTRSNodeIDs.reserve(nbTRSFaces*nbNodesPerFace);
    allTRSNodes.reserve(nbTRSFaces*nbNodesPerFace*m_dimension);
  }
  
  NodeDim nodeDim(CFNULL);
  NodeDim::setSize(m_dimension);
  //  const CFreal eps = std::pow(10.,-1.*m_precision);
  
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
    const CFuint nbTRSFaces = trs->getNbElems();
    // AL: suppose not to have degenerated quadrilaterals (with one duplicated node)
    const CFuint nbNodesPerFace = trs->getNbNodesPerElem();
    for (CFuint iFace = 0; iFace < nbTRSFaces; ++iFace) {
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	totalNbTRSNodes++;
	const CFuint trsNodeID = (*trs)(iFace,iNode);
	trs->setNodeDim(trsNodeID, nodeDim);
	
	// reduce precision to allow for sufficiently high matching tolerance
	// cout << "TRS: " <<  nodeDim << " -> "; 
	// MathTools::MathFunctions::roundUpPrecision(nodeDim, eps);
	///cout << nodeDim << "\n";
	
	
	cf_assert(m_elemMapCoordToNodeID.count(nodeDim) == 1);
	(*trs)(iFace,iNode) = m_elemMapCoordToNodeID.find(nodeDim)->second;
	
	if (m_hasAllSurfFile) {
	  allTRSNodeIDs.push_back((*trs)(iFace,iNode));
	  // store node coordinates with reduced precision
	  for (CFuint iDim = 0; iDim < m_dimension; ++iDim) {
	    allTRSNodes.push_back(nodeDim[iDim]);
	  }
	}
      }
    }
  }
  
  if (m_hasAllSurfFile) {
    cf_assert(totalNbTRSNodes == allTRSNodes.size()/m_dimension);
    cf_assert(allTRSNodes.size() == allTRSNodes.capacity());
    cf_assert(allTRSNodes.size() == allTRSNodeIDs.size()*m_dimension);
  }
  
  CFLog(INFO, "Tecplot2CFmeshConverter::renumberTRSData() => " << counter << "/" << totalNbTRSNodes << " NOT found\n");
  
  // if one single TRS coming from TECPLOT interpolation is available, we now have to 
  // match the TECPLOT-generated boundary points coordinates (.allsurf.plt) with 
  // the COOLFluiD-generated ones (.surf.plt) and assign IDs to the TRS nodes: 
  // to achieve this we rely upon a shortest distance algorithm 
  if (m_hasAllSurfFile) {
    CFuint counter = 0;
    CFuint totalNbTRSNodes = 0;
    const CFuint startTRS = m_nbElemTypes+1;
    const CFuint nbTRSs = m_elementType.size() - startTRS;
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
      const CFuint nbTRSFaces = trs->getNbElems();
      const CFuint nbNodesPerFace = trs->getNbNodesPerElem();
      for (CFuint iFace = 0; iFace < nbTRSFaces; ++iFace) {
	for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	  totalNbTRSNodes++;
	  const CFuint trsNodeID = (*trs)(iFace,iNode);
	  trs->setNodeDim(trsNodeID, nodeDim);
	  
	  // set the nodeID corresponding to the closest boundary point in the current face
	  (*trs)(iFace,iNode) = getClosestNodeID(nodeDim, allTRSNodes, allTRSNodeIDs);
	}
      }
    }
  }
  
  CFLog(VERBOSE, "Tecplot2CFmeshConverter::renumberTRSData() => END/n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFint Tecplot2CFmeshConverter::getClosestNodeID(const NodeDim& nodeDim, 
						const std::vector<CFreal>& allTRSNodes,
						const std::vector<CFuint>& allTRSNodeIDs)
{
  using namespace std;
  using namespace COOLFluiD::MathTools;
  
  vector<CFreal>& allTRSNodesLocal = const_cast<vector<CFreal>&>(allTRSNodes); 
  
  const CFuint nbNodes = allTRSNodeIDs.size();
  NodeDim tmpNode(CFNULL);
  NodeDim::setSize(m_dimension);
  CFint nodeID = -1;
  CFreal distanceSqMin = numeric_limits<CFreal>::max();
  
  for (CFuint i = 0; i < nbNodes; ++i) {
    tmpNode.reset(&allTRSNodesLocal[i*m_dimension]);
    const CFreal distanceSq = MathFunctions::getSquaredDistance(nodeDim, tmpNode);
    if (distanceSq < distanceSqMin) {
      distanceSqMin = distanceSq;
      nodeID = allTRSNodeIDs[i];
    }
  }
  
  cf_assert(nodeID >=0 );
  return nodeID;
}
 
 
//////////////////////////////////////////////////////////////////////////////
    
} // namespace Tecplot2CFmesh

    } // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
