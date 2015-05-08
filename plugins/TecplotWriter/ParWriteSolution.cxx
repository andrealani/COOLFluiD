// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/PE.hh"
#include "Common/MPI/MPIIOFunctions.hh"
#include "Common/CFMap.hh"
#include "Common/OSystem.hh"
#include "Common/BadValueException.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/WriteListMap.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/ParWriteSolution.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

 MethodCommandProvider<ParWriteSolution, TecWriterData, TecplotWriterModule>
 parWriteSolutionProvider("ParWriteSolution");

 //////////////////////////////////////////////////////////////////////////////

 void cmpAndTakeMaxAbs3(CFreal* invec, CFreal* inoutvec, int* len,
			MPI_Datatype* datatype)
 {
   cf_assert(len != CFNULL);
   int size = *len;
   for (int i = 0; i < size; ++i) {
     inoutvec[i] = (fabs(invec[i]) > 0.) ? invec[i] : inoutvec[i];
   }
 }

 //////////////////////////////////////////////////////////////////////////////

 void ParWriteSolution::defineConfigOptions(Config::OptionList& options)
 {
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file."); 
   options.addConfigOption< CFuint >("NbWriters", "Number of writers (and MPI groups)");
   options.addConfigOption< int >("MaxBuffSize", "Maximum buffer size for MPI I/O");
 }

 //////////////////////////////////////////////////////////////////////////////

 ParWriteSolution::ParWriteSolution(const std::string& name) : 
   TecWriterCom(name),
   ParFileWriter(),
   socket_nodes("nodes"),
   socket_nstatesProxy("nstatesProxy"),
   _headerOffset(),
   _oldNbNodesElemsInType(),
   _totalNbNodesInType(),
   _nodesInType(),
   _mapNodeID2NodeIDByEType(),
   _mapGlobal2LocalNodeID()
 {
   addConfigOptionsTo(this);

   _fileFormatStr = "ASCII";
   setParameter("FileFormat",&_fileFormatStr);

   _nbWriters = 1;
   setParameter("NbWriters",&_nbWriters);

   _maxBuffSize = 2147479200; // (CFuint) std::numeric_limits<int>::max();
   setParameter("MaxBuffSize",&_maxBuffSize);
 }

 //////////////////////////////////////////////////////////////////////////////

 ParWriteSolution::~ParWriteSolution() 
 {
   CFAUTOTRACE;
   
   cleanupNodeIDMapping();
 }

 //////////////////////////////////////////////////////////////////////////////

 void ParWriteSolution::execute()
 {
   CFAUTOTRACE;
   
   CFLog(INFO, "Writing solution to " << getMethodData().getFilename().string() << "\n");
   
   if (hasChangedMesh()) {
     buildNodeIDMapping();
   }
   
   if(_fileFormatStr == "ASCII") {    
     // reset to 0 the new file flag
     _isNewFile = 0;
     ofstream* file = CFNULL;
     Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
       Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
     
    if (_isWriterRank) { 
      // if the file has already been processed once, open in I/O mode
      if (_fileList.count(getMethodData().getFilename()) > 0) {
	file = &fhandle->open(getMethodData().getFilename(), ios_base::in | ios_base::out);
      }
      else {
	file = &fhandle->open(getMethodData().getFilename());
	
	// if the file is a new one add it to the file list
	_fileList.insert(getMethodData().getFilename());
	_isNewFile = 1;
      }
    }
    
    int isNewFile = (int)_isNewFile;
    MPI_Bcast(&isNewFile, 1, MPIStructDef::getMPIType(&isNewFile), _ioRank, _comm);
    
    _isNewFile = (isNewFile == 0) ? false : true;
    
    writeToFileStream(getMethodData().getFilename(), file); 
    
    if (_isWriterRank) {
      PE::Group& wg = PE::getGroup("Writers");
      MPI_Barrier(wg.comm);
      fhandle->close();
    }
  }
  else {
    // cf_assert(_fileFormatStr == "BINARY");
    
    // ///@todo change this to use the tecplot library
    // ///this is slow and NOT portable but at least, it takes less space
    // writeToFile("tmp");
    // std::string transformFile = "$TECHOME/bin/preplot tmp " + getMethodData().getFilename().string();
    // CFout << transformFile << "\n";
    
    // OSystem::getInstance().executeCommand(transformFile);
    
    // //     writeToBinaryFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

const std::string ParWriteSolution::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeToBinaryFile()
{
 CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////
      
void ParWriteSolution::writeToFileStream(const boost::filesystem::path& filepath, 
					 std::ofstream* fout)
{
  CFAUTOTRACE;
  
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  
  if (!getMethodData().onlySurface()) {
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    datahandle_output->getDataHandles();
    
    // write the TECPLOT file header
    if (_isNewFile) {
      writeHeader(fout);
    }
    
    // MeshElementType stores GLOBAL element type data
    const vector<MeshElementType>& me =
      MeshDataStack::getActive()->getTotalMeshElementTypes();
    
    std::vector<SafePtr<TopologicalRegionSet> > trsList =
      MeshDataStack::getActive()->getTrsList();
    
    for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs) {
      SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
      
      if ((trs->hasTag("inner")) && (trs->hasTag("cell"))) {
	SafePtr<vector<ElementTypeData> > elementType =
	  MeshDataStack::getActive()->getElementTypeData(trs->getName());
	
	// loop over the element types
	// and create a zone in the tecplot file
	// for each element type
	for (CFuint iType = 0; iType < elementType->size(); ++iType) {
	  ElementTypeData& eType = (*elementType)[iType];
	  const CFuint nbCellsInType  = eType.getNbElems();
	  
	  const bool meshChanged = hasChangedMesh(iType);
	  if (_isNewFile || meshChanged) {
	    vector<MPI_Offset>& headerOffset = _headerOffset[iType];
	    
	    // print zone header: one zone per element type
	    
	    if (_myRank == _ioRank) {
	      headerOffset[0] = fout->tellp();
	      
	      *fout << "ZONE "
		    << "  T= \"ZONE" << iType << " " << eType.getShape() <<"\""
		    << ", N=" << _totalNbNodesInType[iType]
		    << ", E=" << me[iType].elementCount
		    << ", F=FEPOINT"
		    << ", ET=" << MapGeoEnt::identifyGeoEntTecplot
		(eType.getNbNodes(),
		 eType.getGeoOrder(),
		 PhysicalModelStack::getActive()->getDim()) 
		    << ", SOLUTIONTIME=" << subSysStatus->getCurrentTimeDim() << flush;
	      
	      if (getMethodData().getAppendAuxData()) {
		*fout << " AUXDATA TRS=\"" << trs->getName() << "\""
		      << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
		      << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
		      << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
		      << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\"";
	      }
	      *fout << "\n";
	      
	      headerOffset[1] = fout->tellp();
	    }
	    
	    // communicate the current file position to all processor
	    MPI_Bcast(&headerOffset[0], 2, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
	    
	    // backup the total counts for this element type
	    _oldNbNodesElemsInType[iType].first  = _totalNbNodesInType[iType];
	    _oldNbNodesElemsInType[iType].second = me[iType].elementCount;
	  }
	  
	  // print nodal coordinates and stored node variables
	  writeNodeList(fout, iType);
	  
	  if (_isNewFile || meshChanged) {
	    // write element-node connectivity
	    writeElementList(fout, iType, trs);
	  }
	  
	}
      } //end if inner cells
    } //end loop over trs
  } // if only surface
  
  // write boundary surface data
  // writeBoundarySurface();
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::setup() => start\n");
  
  TecWriterCom::setup();
  ParFileWriter::setWriterGroup();
  
  const vector<MeshElementType>& me =
    MeshDataStack::getActive()->getTotalMeshElementTypes();
  const CFuint nbElementTypes = me.size();
  cf_assert(nbElementTypes > 0);
  
  _offset.resize(nbElementTypes);
  
  _headerOffset.resize(nbElementTypes);
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    _headerOffset[i].resize(2, 0);
  }
  
  _oldNbNodesElemsInType.resize(nbElementTypes);
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    _oldNbNodesElemsInType[i].first = 0;
    _oldNbNodesElemsInType[i].second = 0;
  }
  
  _totalNbNodesInType.resize(nbElementTypes, 0);
  
  CFLog(VERBOSE, "ParWriteSolution::setup() => end\n");
}      
      
//////////////////////////////////////////////////////////////////////////////
 
void ParWriteSolution:: cleanupNodeIDMapping()
{  
  CFAUTOTRACE;
  
  for (CFuint i = 0; i < _nodesInType.size(); ++i) {
    vector<CFuint>().swap(_nodesInType[i]);
  }
  
  for (CFuint i = 0; i < _mapNodeID2NodeIDByEType.size(); ++i) {
    delete _mapNodeID2NodeIDByEType[i];
  }
  
  _mapGlobal2LocalNodeID.clear();
}  
 
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution:: buildNodeIDMapping()
{   
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() => start\n");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  // MeshElementType stores GLOBAL element type data
  const vector<MeshElementType>& me = 
    MeshDataStack::getActive()->getTotalMeshElementTypes();
  const CFuint nbElementTypes = me.size();
  
  cleanupNodeIDMapping();
  
  _nodesInType.resize(nbElementTypes);
  _mapNodeID2NodeIDByEType.resize(nbElementTypes);
  _mapGlobal2LocalNodeID.reserve(nodes.size());
  
  vector<SafePtr<TopologicalRegionSet> > trsList = 
    MeshDataStack::getActive()->getTrsList();
  
  const CFuint totalNodeCount =
    MeshDataStack::getActive()->getTotalNodeCount();
  cf_assert(totalNodeCount >= nodes.size());
  vector<bool> foundGlobalID(totalNodeCount);
  
  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
    if ((trs->hasTag("inner")) && (trs->hasTag("cell"))) {
      SafePtr<vector<ElementTypeData> > elementType =
	MeshDataStack::getActive()->getElementTypeData(trs->getName());
      cf_assert(elementType->size() == nbElementTypes);
      
      // AL: check possible inconsistency between MeshElementType and ElementTypeData storages
      for (CFuint iType = 0; iType < elementType->size(); ++iType) {
	ElementTypeData& eType = (*elementType)[iType];
	const CFuint nbCellsInType  = eType.getNbElems();
	
	vector<CFuint>& nodesInType = _nodesInType[iType];
	
	// find which global nodeIDs are used in the elements of this type
	if (nbCellsInType > 0) {
	  const CFuint nbNodesInType  = eType.getNbNodes();
	  nodesInType.reserve(nbCellsInType*nbNodesInType); // this array can be oversized
	  
	  for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) {
	    cf_assert(nbNodesInType == trs->getNbNodesInGeo(iCell));
	    for (CFuint in = 0; in < nbNodesInType; ++in) {
	      nodesInType.push_back(nodes[trs->getNodeID(iCell,in)]->getGlobalID());
	    }
	  }
	  
	  // sort the vector so we can then remove duplicated nodes
	  sort(nodesInType.begin(), nodesInType.end(), std::less<CFuint>());
	  
	  // remove duplicated nodes
	  vector<CFuint>::iterator lastNode = unique(nodesInType.begin(), nodesInType.end());
	  nodesInType.erase(lastNode, nodesInType.end());
	}
	
	// now each process has a unique list of global node IDs for the current element type
	// maximum number of unique nodes across all processors for the current element type
	CFuint maxNbNodesInProc = 0;
	CFuint nbLocalNodesInType = nodesInType.size();
	MPIError::getInstance().check
	  ("MPI_Allreduce", "ParWriteSolution::setup()",
	   MPI_Allreduce(&nbLocalNodesInType, &maxNbNodesInProc, 1, 
			 MPIStructDef::getMPIType(&maxNbNodesInProc), MPI_MAX, _comm));
	
	cf_assert(maxNbNodesInProc >= nbLocalNodesInType);
	
	// temporary buf to store global IDs
	vector<CFuint> buf(maxNbNodesInProc);
	
	// reset the flags for all nodes
	foundGlobalID.assign(foundGlobalID.size(), false);
	
	// use collective broadcast to know all global IDs referencing the current element type
	for (CFuint root = 0; root < _nbProc; ++root) {
	  CFuint bufSize = 0;
	  if (root == _myRank) {
	    for (CFuint i = 0; i < nbLocalNodesInType; ++i) {
	      buf[i] = nodesInType[i];
	    }
	    bufSize = nbLocalNodesInType;
	  } 
	  
	  // communicate the buffer size
	  MPIError::getInstance().check
	    ("MPI_Bcast", "ParWriteSolution::setup()", 
	     MPI_Bcast(&bufSize, 1, MPIStructDef::getMPIType(&bufSize), root, _comm));
	  
	  // the following check is needed for ensuring consistent collective behaviour
	  if (bufSize > 0) {
	    cf_assert(bufSize <= foundGlobalID.size());
	    
	    MPIError::getInstance().check
	      ("MPI_Bcast", "ParWriteSolution::setup()", 
	       MPI_Bcast(&buf[0], bufSize, MPIStructDef::getMPIType(&buf[0]), root, _comm));
	    
	    // flag all global IDs that are found
	    for (CFuint i = 0; i < bufSize; ++i) {
	      foundGlobalID[buf[i]] = true;
	    }
	  }
	}
	
	_mapNodeID2NodeIDByEType[iType] = new CFMap<CFuint, CFuint>(nodesInType.size());
	
	// create mapping between global IDs belonging to the current processor
	// and global IDs reordered by element type (starting from 0) 
	long long int etypeNodeID = -1;
	CFuint countMatching = 0;
	for (CFuint globalNodeID = 0; globalNodeID < totalNodeCount; ++globalNodeID) {
	  if (foundGlobalID[globalNodeID]) {
	    etypeNodeID++; // in reality TECPLOT IDs start from 1 (later on this will have to be considered)
	  }
	  if (std::binary_search(nodesInType.begin(), nodesInType.end(), globalNodeID)) {
	    _mapNodeID2NodeIDByEType[iType]->insert(globalNodeID, etypeNodeID);
	    countMatching++;
	  }
	}
	
	_totalNbNodesInType[iType] = etypeNodeID+1;
	
	cf_assert(countMatching == nodesInType.size());
	
	_mapNodeID2NodeIDByEType[iType]->sortKeys();
      } 
    }
  }
  
  // build the inverse mapping global IDs to local IDs for all nodes 
  for (CFuint i = 0; i < nodes.size(); ++i) {
    _mapGlobal2LocalNodeID.insert(nodes[i]->getGlobalID(), nodes[i]->getLocalID());
  }
  _mapGlobal2LocalNodeID.sortKeys();
  
  CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeBoundarySurface()
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // AL: this is a handle for the nodal states which can be
  // stored as arrays of State*, RealVector* or RealVector
  // but they are always used as arrays of RealVector*
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // we assume that the element-node connectivity
  // is same as element-state connectivity
  /// @todo tecplot writer should not assume element-to-node
  //connectivity to be the same as element-to-state
  //  cf_assert(nodes.size() == nodalStates.getSize());

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;

  const std::vector<std::string>& surfTRS = getMethodData().getSurfaceTRSsToWrite();
  CFuint countTRToWrite = 0;
  std::vector<std::string>::const_iterator itr = surfTRS.begin();

  for(; itr != surfTRS.end(); ++itr) {
    Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);
    const CFuint nbTRs = currTrs->getNbTRs();
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
      SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
      if (tr->getLocalNbGeoEnts() > 0) {
	countTRToWrite++;
      }
    }
  }

  // AL: the file is written only if there is at least one TR with more than 0 nodes
  // else Tecplot cannot handle it and you have manually to skip the file
  if (countTRToWrite > 0) {
    path cfgpath = getMethodData().getFilename();
    path filepath = cfgpath.branch_path() / ( basename(cfgpath) + "-surf" + extension(cfgpath) );

    Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(filepath);

    //  Tecplot Header
    fout << "TITLE      = Boundary data" << "\n";
    fout << "VARIABLES  = ";
    for (CFuint i = 0; i < dim; ++i) {
      fout << " \"x" << i << '\"';
    }

    if (!getMethodData().onlyCoordinates()) {
      SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
      const vector<std::string>& varNames = updateVarSet->getVarNames();
      cf_assert(varNames.size() == nbEqs);

      for (CFuint i = 0 ;  i < nbEqs; ++i) 
      {
	      std::string n = varNames[i];
	      if ( *n.begin()   != '\"' )  n  = '\"' + n;
	      if ( *n.rbegin()  != '\"' )  n += '\"';
	      fout << " " << n;
      }

      if (getMethodData().shouldPrintExtraValues()) {
	vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
	for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
	  fout << " " << extraVarNames[i];
	}
      }
    }

    fout << "\n";

#if 0 // DEBUG
    std::string fileNameTST = getMethodData().getFilename() + ".TEST";
    ofstream tstout(fileNameTST.c_str());
#endif // DEBUG

    // array to store the dimensional states
    RealVector dimState(nbEqs);
    RealVector extraValues; // size will be set in the VarSet
    State tempState;

    std::vector<std::string>::const_iterator itr = surfTRS.begin();
    for(; itr != surfTRS.end(); ++itr) {

      Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);

      const CFuint nbTRs = currTrs->getNbTRs();
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
	SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
	const CFuint nbBFaces = tr->getLocalNbGeoEnts();
	if (nbBFaces > 0) {
	  // get all the unique nodes in the TR
	  vector<CFuint> trNodes;
	  tr->putNodesInTR(trNodes);

	  const CFuint nbAllNodes = trNodes.size();
	  CFMap<CFuint,CFuint> mapNodesID(nbAllNodes);
	  for (CFuint i = 0; i < nbAllNodes; ++i) {
	    mapNodesID.insert(trNodes[i],i+1);
	  }
	  mapNodesID.sortKeys();

	  std::string elemShape;
	  CFuint maxNbNodesInGeo = 2;
	  if (dim == 2) {
	    elemShape = "LINESEG";
	  }
	  else if (dim == 3) {
	    maxNbNodesInGeo = 3;
	    const CFuint nbBFaces = tr->getLocalNbGeoEnts();
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      maxNbNodesInGeo = max(maxNbNodesInGeo, tr->getNbNodesInGeo(iGeo));
	    }

	    // AL: the maximum number of nodes in face is considered:
	    // an extra virtual node is added if TETRA and QUADRILATERAL
	    // are both present
	    elemShape = (maxNbNodesInGeo == 3) ? "TRIANGLE" : "QUADRILATERAL";
	  }

	  // AL: Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
	  if (tr->getLocalNbGeoEnts() > 0) {
	    // print zone header
	    // one zone per TR
	    fout << "ZONE N=" << nbAllNodes
		 << ", T=\"" << currTrs->getName() << ", TR " << iTR << "\""
		 << ", E=" << tr->getLocalNbGeoEnts()
		 << ", F=FEPOINT"
		 << ", ET=" << elemShape
		 << ", SOLUTIONTIME=" << subSysStatus->getCurrentTimeDim()
		 << "\n";
	    
	    vector<CFuint>::const_iterator itr;
	    // print  nodal coordinates and stored nodal variables
	    for (itr = trNodes.begin(); itr != trNodes.end(); ++itr) {
	      // node has to be printed with the right length
	      const CFuint nodeID = *itr;
	      const Node& currNode = *nodes[nodeID];
	      for (CFuint iDim = 0; iDim < dim; ++iDim) {
		fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		fout << currNode[iDim]*refL << " ";
	      }
	      
	      const RealVector& currState = *nodalStates.getState(nodeID);
	      for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
		tempState[ieq] = currState[ieq];
	      }

	      if (!getMethodData().onlyCoordinates()) {
		const CFuint stateID = nodalStates.getStateLocalID(nodeID);
		tempState.setLocalID(stateID);
		SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

		if (getMethodData().shouldPrintExtraValues()) {
		  // dimensionalize the solution
		  updateVarSet->setDimensionalValuesPlusExtraValues
		    (tempState, dimState, extraValues);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << " " << extraValues << "\n";
		}
		else {
		  // set other useful (dimensional) physical quantities
		  updateVarSet->setDimensionalValues(tempState, dimState);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << "\n";
		}
	      }
	      else {
		fout << "\n";
	      }
	    }

	    // #if 0 // DEBUG
	    //      tstout << "FACES"  << std::endl;
	    //         // write  Connectivity
	    //         GeomEntList::const_iterator itg;
	    //         for (itg = bfaces->begin(); itg != bfaces->end(); ++itg) {
	    //           const vector<Node*>& nodesInFace = *(*itg)->getNodes();
	    //           vector<Node*>::const_iterator itr;
	    //           for (itr = nodesInFace.begin(); itr != nodesInFace.end();
	    //            ++itr) {
	    //             tstout << (*itr)->getGlobalID() << "-" << (*itr)->getLocalID()
	    //                 << " : ";
	    //           }
	    //           tstout << std::endl;
	    //         }
	    //      tstout << "CONNECTIVITY"  << std::endl;
	    // #endif

	    // write  Connectivity
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      const CFuint nbNodesInGeo = tr->getNbNodesInGeo(iGeo);
	      for (CFuint in = 0; in < nbNodesInGeo; ++in) {
#if 0 // DEBUG
		tstout << tr->getNodeID(iGeo,in); tstout.flush();
		tstout << " : " << mapNodesID.find(tr->getNodeID(iGeo,in)) << endl;
#endif
		fout << mapNodesID.find(tr->getNodeID(iGeo,in)) << " ";
	      }
	      // if the number of face nodes is less than the
	      // maximum number of nodes in a face of this TR
	      // add an extra dummy node equal to the last "real" one
	      if (nbNodesInGeo < maxNbNodesInGeo) {
		cf_assert(maxNbNodesInGeo == nbNodesInGeo + 1);
		fout << mapNodesID.find(tr->getNodeID(iGeo, (nbNodesInGeo-1))) << " ";
	      }
	      fout << "\n";
	      fout.flush();
	    }
	  }
	}
      }
    }

    fout.close();
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ParWriteSolution::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_nstatesProxy);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
      
void ParWriteSolution::writeHeader(MPI_File* fh)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() start\n");
  
  if (_myRank == _ioRank) {
    
    //  Tecplot Header
    // maximum string size to write includes 30 characters 
    MPIIOFunctions::writeKeyValue<char>(fh, "TITLE = ");
    MPIIOFunctions::writeKeyValue<char>(fh, "Unstructured grid data");
    MPIIOFunctions::writeKeyValue<char>(fh, "\nVARIABLES = ");
    
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < dim; ++i) {
      MPIIOFunctions::writeKeyValue<char>(fh, " \"x" + StringOps::to_str(i) + '\"');
    }
    
    SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
    const vector<std::string>& varNames = updateVarSet->getVarNames();
    cf_assert(varNames.size() == nbEqs);
    
    for (CFuint i = 0 ;  i < nbEqs; ++i)  {
      std::string n = varNames[i];
      if ( *n.begin()  != '\"' )  n  = '\"' + n;
      if ( *n.rbegin() != '\"' )  n += '\"';
      MPIIOFunctions::writeKeyValue<char>(fh, " " + n);
    }
    
    if (getMethodData().shouldPrintExtraValues()) {
      vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
      for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
    	MPIIOFunctions::writeKeyValue<char>(fh, " " +  extraVarNames[i]);
      }
    }
    
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    vector<string> dhVarNames = datahandle_output->getVarNames();
    for (CFuint i = 0 ;  i < dhVarNames.size(); ++i) {
      MPIIOFunctions::writeKeyValue<char>(fh, " " + dhVarNames[i]);
    }
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() end\n");
}
  
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeHeader(std::ofstream* fout) 
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() start\n");
  
  if (_myRank == _ioRank) {
    //  Tecplot Header
    // maximum string size to write includes 30 characters 
    *fout << "TITLE = Unstructured grid data";
    *fout << "\nVARIABLES = ";
    
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < dim; ++i) {
      *fout << " \"x" << StringOps::to_str(i) << '\"';
    }
    
    SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
    const vector<std::string>& varNames = updateVarSet->getVarNames();
    cf_assert(varNames.size() == nbEqs);
    
    for (CFuint i = 0 ;  i < nbEqs; ++i)  {
      std::string n = varNames[i];
      if ( *n.begin()  != '\"' )  n  = '\"' + n;
      if ( *n.rbegin() != '\"' )  n += '\"';
      *fout << " " << n;
    }
    
    if (getMethodData().shouldPrintExtraValues()) {
      vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
      for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
	*fout << " " <<  extraVarNames[i];
      }
    }
    
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    vector<string> dhVarNames = datahandle_output->getVarNames();
    for (CFuint i = 0 ;  i < dhVarNames.size(); ++i) {
      *fout << " " << dhVarNames[i];
    }
    *fout << "\n";
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() end\n");
}
  
//////////////////////////////////////////////////////////////////////////////
  
void ParWriteSolution::writeNodeList(ofstream* fout, const CFuint iType)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeNodeList() => start\n");
    
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  
  CFuint nodesStride = dim + nbEqs;
  if (getMethodData().shouldPrintExtraValues()) {
    nodesStride += updateVarSet->getExtraVarNames().size();
  }
  nodesStride += datahandle_output->getVarNames().size();
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);
  
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // this isa sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = 1; // just nodes
  const CFuint nbLocalElements = _nodesInType[iType].size();
  const CFuint totNbNodes = _totalNbNodesInType[iType];
  RealVector dimState(nbEqs);
  RealVector extraValues; // size will be set in the VarSet
  State tempState;
  
  PE::Group& wg = PE::getGroup("Writers");
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);
  
  // fill in the writer ist
  CFuint totalToSend = 0;
  elementList.fill(totNbNodes, nodesStride, totalToSend);
  
  const CFuint wordFormatSize = 22;
  
  vector<MPI_Offset> wOffset(_nbWriters, _headerOffset[iType][1]); 
  // set the offsets for the nodes of this element type
  _offset[iType].nodes.first  = _headerOffset[iType][1];
  _offset[iType].nodes.second = _offset[iType].nodes.first + totalToSend*wordFormatSize;
  
  CFLog(VERBOSE, "ParWriteSolution::writeNodeList() => offsets = [" 
	<<  _offset[iType].nodes.first << ", " << _offset[iType].nodes.second << "]\n");
  
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = elementList.getMaxElemSize();
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  const vector<CFuint>& nodesInType = _nodesInType[iType];
  for (CFuint iElem = 0; iElem < nbLocalElements; ++iElem) {
    const CFuint globalTypeID = _mapNodeID2NodeIDByEType[iType]->find(nodesInType[iElem]);
    elementList.insertElemLocalID(iElem, globalTypeID, 0);
  }
  elementList.endElemInsertion(_myRank);
  
  // buffer data to send
  vector<CFreal> sendElements(maxElemSendSize, 0);
  vector<CFreal> elementToPrint(maxElemSendSize, 0);
  
  CFint wRank = -1; 
  CFuint wSendSize = 0;
  CFuint rangeID = 0;
  CFuint countElem = 0;
  for (CFuint is = 0; is < nSend; ++is, ++rangeID) {
    bool isRangeFound = false;
    WriteListMap::List elist = elementList.find(rangeID, isRangeFound);
    
    if (isRangeFound) {
      CFuint eSize = 0;
      for (WriteListMap::ListIterator it = elist.first; it != elist.second; ++it, ++eSize) {
	const CFuint localElemID = it->second;
	const CFuint globalElemID = _mapNodeID2NodeIDByEType[iType]->find(nodesInType[localElemID]);
	const CFuint sendElemID = globalElemID - countElem;
	const CFuint nodeID = _mapGlobal2LocalNodeID.find(nodesInType[localElemID]);
	
	// this fix has to be added EVERYWHERE when writing states in parallel
	if (nodes[nodeID]->isParUpdatable()) {
	  CFuint isend = sendElemID*nodesStride;
	  for (CFuint in = 0; in < dim; ++in, ++isend) {
	    cf_assert(isend < sendElements.size());
	    cf_assert(nodeID < nodes.size());
	    sendElements[isend] = (*nodes[nodeID])[in]*refL;
	  }
	  
	  const RealVector& currState = *nodalStates.getState(nodeID);
	  for (CFuint in = 0; in < nbEqs; ++in, ++isend) {
	    cf_assert(isend < sendElements.size());
	    sendElements[isend] = currState[in];
	  }
	  
	  const CFuint stateID = nodalStates.getStateLocalID(nodeID);
	  tempState.setLocalID(stateID);
	  // the node is set  in the temporary state
	  tempState.setSpaceCoordinates(nodes[nodeID]);
	  
	  if (getMethodData().shouldPrintExtraValues()) {
	    // dimensionalize the solution
	    updateVarSet->setDimensionalValuesPlusExtraValues
	      (tempState, dimState, extraValues);
	    for (CFuint in = 0; in < extraValues.size(); ++in, ++isend) {
	      cf_assert(isend < sendElements.size());
	      sendElements[isend] = extraValues[in];
	    }
	  }
	  
	  datahandle_output->fillStateData(&sendElements[0], stateID, isend);
	}
      }
      
      cf_assert(eSize*nodesStride <= elementList.getSendDataSize(rangeID));
    }
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" sendElements  = ", &sendElements, nodesStride) << "\n");
    
    const CFuint sendSize = elementList.getSendDataSize(rangeID);
    cf_assert(sendSize <= sendElements.size());
    cf_assert(sendSize <= elementToPrint.size());
    
    // if the rank corresponds to a writing process, record the size to send for this writer
    if (_isWriterRank && wg.globalRanks[is] == _myRank) {
      wSendSize = sendSize; // this should be the total sendsize in the range
      wRank = is;
    }
    
    MPI_Op myMpiOp;
    MPI_Op_create((MPI_User_function *)cmpAndTakeMaxAbs3, 1, &myMpiOp);
    MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
	       MPIStructDef::getMPIType(&sendElements[0]), myMpiOp, wg.globalRanks[is], _comm);
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" elementToPrint  = ", &elementToPrint, nodesStride) << "\n");
    
    if (is == 0) {
      CFLog(VERBOSE, "[0] => wOffset = " << wOffset[0] << ", sendSize = " << sendSize << "\n");
    }
    
    // the offsets for all writers with send ID > current must be incremented  
    for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
      cf_assert(sendSize > 0);
      wOffset[iw] += sendSize*wordFormatSize;
      CFLog(VERBOSE, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
    }
    
    // reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      sendElements[i] = 0;
    }
    
    // update the count element for the current element type
    countElem += elementList.getSendDataSize(rangeID)/nodesStride;
  }
  
  if (_isWriterRank) { 
    CFLog(DEBUG_MIN, "ParWriteSolution::writeNodeList() => P[" << _myRank << "] => offset = " << wOffset[wRank] << "\n");
    // each writer can now concurrently write all the collected data (related to one element type)
    cf_assert(wRank >= 0);
    
    // point to the corresponding writing location
    fout->seekp(wOffset[wRank]);
    
    CFuint countN = 0;
    CFLog(VERBOSE, "wSendSize = " << wSendSize << ", nodesStride = " << nodesStride << "\n");
    for (CFuint i = 0; i < wSendSize; ++i) {
      const CFreal coeff = (countN < dim) ? refL : 1.;
      // this format corresponds to a line of 22*nodesStride bytes  
      fout->precision(14);
      fout->setf(ios::scientific,ios::floatfield); 
      fout->setf(ios::showpos);
      if ((i+1)%nodesStride > 0) {
	*fout << elementToPrint[i]*coeff << " ";
	countN++;
      }
      else {
	*fout << elementToPrint[i]*coeff << "\n";
	countN = 0; // reset countN to 0
      }
    }
    
    MPI_Offset lastpos = fout->tellp();
    MPI_Offset maxpos  = 0;
    MPI_Allreduce(&lastpos, &maxpos, 1, MPIStructDef::getMPIOffsetType(), MPI_MAX, wg.comm);
    
    cf_assert(_offset[iType].nodes.second == maxpos);
  }
  
  //reset the all sendElement list to 0
  for (CFuint i = 0; i < maxElemSendSize; ++i) {
    elementToPrint[i] = 0;
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeNodeList() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeElementList(ofstream* fout,
					const CFuint iType,
					SafePtr<TopologicalRegionSet> elements)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeElementList() => start\n");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  // MeshElementType stores GLOBAL element type data
  const vector<MeshElementType>& me =
    MeshDataStack::getActive()->getTotalMeshElementTypes();
  
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();
  const ElementTypeData& eType = (*elementType)[iType];
  
  cf_assert(elementType->size() == me.size());
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = 1; // one element type at a time
  const CFuint nbLocalElements = eType.getNbElems();
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);
  
  // store global info about the global ID ranges for sending
  // and the size of each send
  // each send will involve one writing process which wil collect all data 
  const CFuint nbElementsInType = me[iType].elementCount;
  const CFuint nbNodesInType    = me[iType].elementNodes;
  const CFuint totalSize = nbElementsInType*nbNodesInType;
  
  // fill in the writer list
  CFuint totalToSend = 0;
  elementList.fill(nbElementsInType, nbNodesInType, totalToSend);
  cf_assert(totalToSend == totalSize);
  
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = elementList.getMaxElemSize();
  const CFuint wordFormatSize = 15;
  
  // wOffset is initialized with current offset
  MPI_Offset offset = _offset[iType].nodes.second;
  vector<MPI_Offset> wOffset(_nbWriters, offset); 
  
  // start element list offset (current position)
  _offset[iType].elems.first = offset;
  // end element list offset
  _offset[iType].elems.second = _offset[iType].elems.first + totalSize*wordFormatSize;
  
  CFLog(VERBOSE, "ParWriteSolution::writeElementList() => offsets = [" 
	<<  _offset[iType].elems.first << ", " << _offset[iType].elems.second << "]\n");
  
  SafePtr< vector<CFuint> > globalElementIDs = 
    MeshDataStack::getActive()->getGlobalElementIDs();
  cf_assert(globalElementIDs->size() >= nbLocalElements);
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  CFuint elemID = eType.getStartIdx();
  const CFuint nbLocalElementsInType = eType.getNbElems();
  for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem, ++elemID) {
    elementList.insertElemLocalID(elemID, (*globalElementIDs)[elemID], 0);
  }
  elementList.endElemInsertion(_myRank);
  
  cf_assert(elemID == eType.getEndIdx());
  CFLog(VERBOSE,_myRank << " maxElemSendSize = " << maxElemSendSize << "\n");
  
  // buffer data to send
  vector<CFuint> sendElements(maxElemSendSize, 0);
  vector<CFuint> elementToPrint(maxElemSendSize, 0);
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  PE::Group& wg = PE::getGroup("Writers");
  CFint wRank = -1; 
  CFuint rangeID = 0;
  CFuint wSendSize = 0;
  CFuint countElem = 0;
  for (CFuint is = 0; is < nSend; ++is, ++rangeID) {
    bool isRangeFound = false;
    WriteListMap::List elist = elementList.find(rangeID, isRangeFound);
    
    if (isRangeFound) {
      CFuint eSize = 0;
      for (WriteListMap::ListIterator it = elist.first; it != elist.second; ++it, ++eSize) {
	const CFuint localElemID = it->second;
	const CFuint globalElemID = (*globalElementIDs)[localElemID];
	const CFuint sendElemID = globalElemID - countElem;
	
	CFuint isend = sendElemID*nbNodesInType;
	for (CFuint in = 0; in < nbNodesInType; ++in, ++isend) {
	  const CFuint localNodeID = elements->getNodeID(localElemID, in);
	  if (isend >= sendElements.size()) {
	    CFLogInfo(_myRank << " nbNodesInType = " << nbNodesInType
		      << " node isend = " << isend << " , size = " << sendElements.size() << "\n");
	    cf_assert(isend < sendElements.size());
	  }
	  const CFuint globalNodeID = nodes[localNodeID]->getGlobalID();
	  sendElements[isend] = _mapNodeID2NodeIDByEType[iType]->find(globalNodeID);
	}
      }
      
      cf_assert(eSize*nbNodesInType <= elementList.getSendDataSize(rangeID));
    }
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		  (" sendElements  = ", &sendElements, nbNodesInType) << "\n");
    
    // size of the data that, in total, those processes will send to the aggregator (writer) with ID=is
    const CFuint sendSize = elementList.getSendDataSize(rangeID);
    cf_assert(sendSize <= sendElements.size());
    cf_assert(sendSize <= elementToPrint.size());
    
    // if the rank corresponds to a writing process, record the size to send for this writer
    if (_isWriterRank && wg.globalRanks[is] == _myRank) {
      wSendSize = sendSize; // this should be the total sendsize in the range
      wRank = is;
    }
    
    // for each send, accumulate data to the corresponding writing process
    // this might allocate too much memory on master processes with OPENMPI
    MPI_Reduce(&sendElements[0], &elementToPrint[0], (int)sendSize,
	       MPIStructDef::getMPIType(&sendElements[0]), MPI_MAX, wg.globalRanks[is], _comm);
    
    // the offsets for all writers with send ID > current must be incremented  
    for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
      cf_assert(sendSize > 0);
      wOffset[iw] += sendSize*wordFormatSize;
      CFLog(DEBUG_MIN, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
    }
    
    //reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      sendElements[i] = 0;
    }
    
    // update the count element for the current element type
    countElem += elementList.getSendDataSize(rangeID)/nbNodesInType;
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		  (" elementToPrint  = ", &elementToPrint, nbNodesInType) << "\n");
    
  } // end sending loop
  
  if (_isWriterRank) { 
    CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeElementList() => P[" << _myRank 
	  << "] => offset = " << wOffset[wRank] << "\n");
    
    // each writer can now concurrently write all the collected data (related to one element type)
    cf_assert(wRank >= 0);
    
    // point to the corresponding writing location
    fout->seekp(wOffset[wRank]);
    
    CFuint countN = 0;
    const CFuint nbElementsToWrite = wSendSize/nbNodesInType;
    for (CFuint i = 0; i < nbElementsToWrite; ++i) {
      writeElementConn(*fout, &elementToPrint[i*nbNodesInType],
		       nbNodesInType, eType.getGeoOrder(), dim);
      *fout << "\n";
    }
    
    MPI_Offset lastpos = fout->tellp();
    MPI_Offset maxpos  = 0;
    MPI_Allreduce(&lastpos, &maxpos, 1, MPIStructDef::getMPIOffsetType(), MPI_MAX, wg.comm);
    
    cf_assert(_offset[iType].elems.second == maxpos);
  }
  
  // reset all the elements to print to 0
  for (CFuint i = 0; i < maxElemSendSize; ++i) {
    elementToPrint[i] = 0;
  }
    
  CFLog(VERBOSE, "ParWriteSolution::writeElementList() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeElementConn(ofstream& file,
					CFuint* nodeIDs,
					const CFuint nbNodes,
					const CFuint geoOrder,
					const CFuint dim)
{
  /// @todo only 1st order geometry
  cf_assert(geoOrder == CFPolyOrder::ORDER1);
  
  switch(dim) {
  case DIM_2D:
    switch(nbNodes) {
      
    case 3:  // TRIANGLE
      
      file << std::noshowpos 
	   << setw(14) << nodeIDs[0]+1 << " "
	   << setw(14) << nodeIDs[1]+1 << " "
	   << setw(14) << nodeIDs[2]+1;
      
      break;
      
    case 4: // QUADRILATERAL
      
      file << std::noshowpos 
	   << setw(14) << nodeIDs[0]+1 << " "
	   << setw(14) << nodeIDs[1]+1 << " "
	   << setw(14) << nodeIDs[2]+1 << " "
	   << setw(14) << nodeIDs[3]+1;
      
      break;
    default:
      std::string msg = std::string("Wrong number of nodes in 2D: ") +
	Common::StringOps::to_str(nbNodes);
      throw BadValueException(FromHere(),msg);
    }
    break;
  case DIM_3D:
    switch(nbNodes) {
      
    case 4: // TETRAHEDRON
      
      file << std::noshowpos 
	   << setw(14) << nodeIDs[0]+1 << " "
	   << setw(14) << nodeIDs[1]+1 << " "
	   << setw(14) << nodeIDs[2]+1 << " "
	   << setw(14) << nodeIDs[3]+1;
      
      break;

    case 5: // BRICK with nodes coalesced 5,6,7->4

      file << std::noshowpos 
	   << setw(14) << nodeIDs[0]+1 << " "
	   << setw(14) << nodeIDs[1]+1 << " "
	   << setw(14) << nodeIDs[2]+1 << " "
	   << setw(14) << nodeIDs[3]+1 << " "
	   << setw(14) << nodeIDs[4]+1 << " "
	   << setw(14) << nodeIDs[4]+1 << " "
	   << setw(14) << nodeIDs[4]+1 << " "
	   << setw(14) << nodeIDs[4]+1;
      
      break;

    case 6: // BRICK with nodes 2->3 and 6->7 coalesced
      
      file << std::noshowpos 
	   << setw(14) << nodeIDs[0]+1 << " "
	   << setw(14) << nodeIDs[1]+1 << " "
	   << setw(14) << nodeIDs[2]+1 << " "
	   << setw(14) << nodeIDs[2]+1 << " "
	   << setw(14) << nodeIDs[3]+1 << " "
	   << setw(14) << nodeIDs[4]+1 << " "
	   << setw(14) << nodeIDs[5]+1 << " "
	   << setw(14) << nodeIDs[5]+1;
      
      break;
      
    case 8: // BRICK
      
      file << std::noshowpos 
	   << setw(14) << nodeIDs[0]+1 << " "
	   << setw(14) << nodeIDs[1]+1 << " "
	   << setw(14) << nodeIDs[2]+1 << " "
	   << setw(14) << nodeIDs[3]+1 << " "
	   << setw(14) << nodeIDs[4]+1 << " "
	   << setw(14) << nodeIDs[5]+1 << " "
	   << setw(14) << nodeIDs[6]+1 << " "
	   << setw(14) << nodeIDs[7]+1;
      
      break;

    default:
      std::string msg = std::string("Wrong number of nodes in 3D: ") +
                   Common::StringOps::to_str(nbNodes);
      throw BadValueException(FromHere(),msg);
    }
    break;
  default:
    std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
      Common::StringOps::to_str(dim);
    throw BadValueException(FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////
      
bool ParWriteSolution::hasChangedMesh() const
{
  CFAUTOTRACE;
  
  // use signals monitoring for this in the future
  for (CFuint iType = 0; iType < _totalNbNodesInType.size(); ++iType) {
    if (hasChangedMesh(iType)) {
      CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh() => true\n");
      return true;
    }
  }
  CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh() => false\n");
  return false;
}
      
//////////////////////////////////////////////////////////////////////////////

bool ParWriteSolution::hasChangedMesh(const CFuint iType) const
{
  CFAUTOTRACE;
  
  // use signals monitoring for this in the future
  cf_assert(iType < _totalNbNodesInType.size());
  cf_assert(iType < _oldNbNodesElemsInType.size());
  if (_totalNbNodesInType[iType] != _oldNbNodesElemsInType[iType].first) {
    CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh (iType = " << iType << ") => true\n");
    return true;
  }
  
  if (MeshDataStack::getActive()->getTotalMeshElementTypes()[iType].elementCount 
      != _oldNbNodesElemsInType[iType].second) {
    CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh (iType = " << iType << ") => true\n");
    return true;
  }
  
  CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh (iType = " << iType << ") => false\n");
  return false;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter
  
} // namespace COOLFluiD
