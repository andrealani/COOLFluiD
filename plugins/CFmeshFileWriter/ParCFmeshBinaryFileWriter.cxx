// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>
#include <numeric>

#include "CFmeshFileWriter/ParCFmeshBinaryFileWriter.hh"
#include "Framework/ElementTypeData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/ElementTypeData.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Common/CFMultiMap.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/MPI/MPIIOFunctions.hh"

#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/WriteListMap.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

    namespace CFmeshFileWriter {


//////////////////////////////////////////////////////////////////////////////

void cmpAndTakeMaxAbs2(CFreal* invec, CFreal* inoutvec, int* len,
		       MPI_Datatype* datatype)
{
  cf_assert(len != CFNULL);
  int size = *len;
  for (int i = 0; i < size; ++i) {
    inoutvec[i] = (fabs(invec[i]) > 0.) ? invec[i] : inoutvec[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

ParCFmeshBinaryFileWriter::ParCFmeshBinaryFileWriter() :
  ParFileWriter(), 
  ConfigObject("ParCFmeshBinaryFileWriter"),
  _writeData()
{ 
  addConfigOptionsTo(this);
  
  _nbWriters = 1;
  setParameter("NbWriters",&_nbWriters);
  
  _maxBuffSize = 2147479200; // (CFuint) std::numeric_limits<int>::max();
  setParameter("MaxBuffSize",&_maxBuffSize);
}
      
//////////////////////////////////////////////////////////////////////////////

ParCFmeshBinaryFileWriter::~ParCFmeshBinaryFileWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("NbWriters", "Number of writers (and MPI groups)");
  options.addConfigOption< int >("MaxBuffSize", "Maximum buffer size for MPI I/O");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::setup()
{
  ParFileWriter::setWriterGroup();
  _offset.resize(1);
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeToFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeToFile() => start\n");
  
  /* MPI_Offset offset = sizeof(double)*pstart;
     MPI_File_set_view( fh0, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL ) ; 
     MPI_File_set_atomicity( fh0, TRUE ) ; 
     MPI_File_write_at(fh0, 0, a, 10, MPI_INT, &status) ; 
     MPI_File_seek(file, offset, MPI_SEEK_SET); */
  
  // AL: the writer processor with less elements is chosen as the master IO node
  // this can avoid memory "explosion" on master node in cases for which 
  // the partitioning is not well balanced
  // make sure that only writer ranks contribute to the selection by assigning 
  // huge number of elements to the oher processors 
  CFuint nbLocalElements = (_isWriterRank) ? getWriteData().getNbElements() : std::numeric_limits<CFuint>::max();
  CFuint minNumberElements = 0;
  MPI_Allreduce(&nbLocalElements, &minNumberElements, 1, MPIStructDef::getMPIType(&nbLocalElements), MPI_MIN, _comm);
  CFuint rank = (minNumberElements == nbLocalElements)  ? _myRank : 0;
  // IO rank is maximum rank whose corresponding process has minimum number of elements
  MPI_Allreduce(&rank, &_ioRank, 1, MPIStructDef::getMPIType(&rank), MPI_MAX, _comm);    
  CFLog(INFO, "ParCFmeshBinaryFileWriter::writeToFile() => IO rank is " << _ioRank << "\n");
  
  // if the file has already been processed once, open in I/O mode
  if (_fileList.count(filepath) == 0) {
    // add the new file to the list
    _fileList.insert(filepath);
    _isNewFile = true;
  }
  
  if (_myRank == _ioRank && (!_isWriterRank)) {
    CFLog(ERROR, "ERROR: ParCFmeshBinaryFileWriter::writeToFile() => IO rank is not a writer rank!\n"); abort();
  }
  
  writeToFileStream(filepath); 
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeToFile() => end\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void ParCFmeshBinaryFileWriter::writeToFileStream
(const boost::filesystem::path& filepath)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeToFileStream() start\n");
   
  // last group is writers group
  PE::Group& wg = PE::getGroup("Writers");
  
  char* fileName = const_cast<char*>(filepath.string().c_str()); 
  
  CFLog(VERBOSE, "fileName = " << filepath.string() << "\n");
  CFLog(VERBOSE, "wg.globalRanks.size() = " << wg.globalRanks.size() << "\n");
  CFLog(VERBOSE, "wg.groupRanks.size() = " << wg.groupRanks.size() << "\n");  
  // all writers open the file for the second or more time
  if (_isWriterRank) {
    MPI_File_open(wg.comm, fileName, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &_fh); 
  }
  
  if (_isNewFile) {
    writeVersionStamp(&_fh);
    
    // global counts
    writeGlobalCounts(&_fh);
    
    // extra vars info
    writeExtraVarsInfo(&_fh);
    
    if (_isWriterRank) {
      MPI_Barrier(wg.comm);
    }
    
    // elements info
    writeElements(&_fh);
    
    // TRS data
    writeTrsData(&_fh);
    
    if (_isWriterRank) {
      _mapFileToStartNodeList[filepath] = _offset[0].TRS.back().second;
    }
  }
  else {
    // what follows works but can leave some lines at the end of the file
    // due to previous writing ...
    cf_assert(_isNewFile == 0);
    
    MPI_Offset position = 0;
    if (_isWriterRank) {
      position = _mapFileToStartNodeList.find(filepath)->second;
      MPI_File_seek(_fh, position, MPI_SEEK_SET);
    }
    
    // communicate to all the processors about the position in the current file
    MPI_Bcast(&position, 1, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
  }
  
  // extra vars that are not state or node related
  writeExtraVars(&_fh);
  
  // write the node list
  writeNodeList(&_fh);
  
  // write the state list
  writeStateList(&_fh);
  
  // terminate the file
  writeEndFile(&_fh);
  
  if (_isWriterRank) {
    MPI_File_close(&_fh);
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeFile() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeVersionStamp(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeDimension() start\n");
  
  if (_myRank == _ioRank) {
    MPIIOFunctions::writeKeyValue<char>(fh, "!COOLFLUID_VERSION ");
    MPIIOFunctions::writeKeyValue<char>(fh, CFEnv::getInstance().getCFVersion());
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!COOLFLUID_SVNVERSION ");
    MPIIOFunctions::writeKeyValue<char>(fh, CFEnv::getInstance().getSvnVersion());
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!CFMESH_FORMAT_VERSION ");
    MPIIOFunctions::writeKeyValue<char>(fh, "1.3");
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeDimension() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeGlobalCounts(MPI_File* fh)
{
 CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeDimension() start\n");
 
 if (_myRank  == _ioRank) {
   MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_DIM ", false, PhysicalModelStack::getActive()->getDim());
   MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_EQ ", false, PhysicalModelStack::getActive()->getNbEq());
   
   MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_NODES ", false, MeshDataStack::getActive()->getTotalNodeCount());
   CFuint nuNodes = getWriteData().getNbNonUpdatableNodes();
   MPI_File_write(*fh, &nuNodes, 1, MPIStructDef::getMPIType(&nuNodes), &_status);
   
   MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_STATES ", false, MeshDataStack::getActive()->getTotalStateCount());
   CFuint nuStates = getWriteData().getNbNonUpdatableStates();
   MPI_File_write(*fh, &nuStates, 1, MPIStructDef::getMPIType(&nuStates), &_status);
   
   const std::vector<CFuint>& tElem = MeshDataStack::getActive()->getTotalElements();
   MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_ELEM ", false, (CFuint)std::accumulate(tElem.begin(), tElem.end(), 0));
 }
 
 CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeDimension() end\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void ParCFmeshBinaryFileWriter::writeExtraVarsInfo(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeExtraVarsInfo() start\n");

 if (_myRank  == _ioRank) {
   if(getWriteData().storePastStates()) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!STORE_PASTSTATES ", false, getWriteData().storePastStates());
   }
   
   if(getWriteData().storePastNodes()) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!STORE_PASTNODES ", false, getWriteData().storePastNodes());
   }

   if(getWriteData().storeInterStates()) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!STORE_INTERSTATES ", false, getWriteData().storeInterStates());
   }
   
   if(getWriteData().storeInterNodes()) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!STORE_INTERNODES ", false, getWriteData().storeInterNodes());
   }
   
   const CFuint nbExtraNodalVars = getWriteData().getNbExtraNodalVars();
   if(nbExtraNodalVars > 0) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_EXTRA_NVARS ", false, nbExtraNodalVars);
     
     MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_NVARS_NAMES ");
     for(CFuint iVar = 0; iVar < nbExtraNodalVars; iVar++) {
       MPIIOFunctions::writeKeyValue<char>(fh, (*(getWriteData().getExtraNodalVarNames()))[iVar] + " ");
     }
     MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_NVARS_STRIDES ");
     MPI_File_write(*fh, &(*(getWriteData().getExtraNodalVarStrides()))[0], 
		    (int)nbExtraNodalVars, 
		    MPIStructDef::getMPIType(&(*(getWriteData().getExtraNodalVarStrides()))[0]), &_status); 
   }
   
   const CFuint nbExtraStateVars = getWriteData().getNbExtraStateVars();
   if (nbExtraStateVars > 0) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_EXTRA_SVARS ", false, nbExtraStateVars);
     
     MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_SVARS_NAMES ");
     for(CFuint iVar = 0; iVar < nbExtraStateVars; iVar++) {
       MPIIOFunctions::writeKeyValue<char>(fh, (*(getWriteData().getExtraStateVarNames()))[iVar] + " ");
     }
     MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_SVARS_STRIDES ");
     MPI_File_write(*fh, &(*(getWriteData().getExtraStateVarStrides()))[0], 
		    (int)nbExtraStateVars, 
		    MPIStructDef::getMPIType(&(*(getWriteData().getExtraStateVarStrides()))[0]), &_status);
   }
   
   const CFuint nbExtraVars = getWriteData().getNbExtraVars();
   if (nbExtraVars > 0) {
     MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_EXTRA_VARS ", false, nbExtraVars);
     
     MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_VARS_NAMES ");
     for(CFuint iVar = 0; iVar < nbExtraVars; iVar++) {
       MPIIOFunctions::writeKeyValue<char>(fh, (*(getWriteData().getExtraVarNames()))[iVar] + " ");
     }
     MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_VARS_STRIDES ");
     MPI_File_write(*fh, &(*(getWriteData().getExtraVarStrides()))[0], 
		    (int)nbExtraVars, 
		    MPIStructDef::getMPIType(&(*(getWriteData().getExtraVarStrides()))[0]), &_status);
   }
 }
 
 CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeExtraVarsInfo() end\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParCFmeshBinaryFileWriter::writeExtraVars(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeExtraVars() start\n");
  
  if (_myRank == _ioRank) {
    getWriteData().prepareExtraVars();
    const vector<CFuint>& extraVarsStrides = *getWriteData().getExtraVarStrides();
    const CFuint totalNbExtraVars = std::accumulate(extraVarsStrides.begin(),extraVarsStrides.end(), 0);
    RealVector& extraValues = getWriteData().getExtraValues();
    cf_assert(extraValues.size() == totalNbExtraVars);
    
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!EXTRA_VARS ");
    if (extraValues.size() > 0) {
      MPI_File_write(*fh, &extraValues[0], (int)totalNbExtraVars, 
		     MPIStructDef::getMPIType(&extraValues[0]), &_status);
    }
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeExtraVars() end\n");
}  

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeElements(MPI_File* fh)
{
  CFLog(VERBOSE,"ParCFmeshBinaryFileWriter::writeElements() start\n");
    
  if (_myRank  == _ioRank) {
    MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_ELEM_TYPES ", 
			  false, MeshDataStack::getActive()->getTotalElements().size());
    MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!GEOM_POLYORDER ", 
			  false, static_cast<CFuint> (getWriteData().getGeometricPolyOrder()));
    MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!SOL_POLYORDER ", 
			  false, static_cast<CFuint> (getWriteData().getSolutionPolyOrder()));
    
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!ELEM_TYPES ");
    vector<MeshElementType>& me = MeshDataStack::getActive()->getTotalMeshElementTypes();
    const CFuint nbElementTypes = me.size();
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      MPIIOFunctions::writeKeyValue<char>(fh, me[i].elementName + " ");
    }
    
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!NB_ELEM_PER_TYPE ");
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      MPI_File_write(*fh, &me[i].elementCount, 1, 
		     MPIStructDef::getMPIType(&me[i].elementCount), &_status); 
    }
    
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!NB_NODES_PER_TYPE ");
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      MPI_File_write(*fh, &me[i].elementNodes, 1, 
		     MPIStructDef::getMPIType(&me[i].elementNodes), &_status); 
    }
    
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!NB_STATES_PER_TYPE ");
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      MPI_File_write(*fh, &me[i].elementStates, 1,
		     MPIStructDef::getMPIType(&me[i].elementStates), &_status); 
    }
  }
  
  // write the list of elements
  writeElementList(fh);

  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeElements() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeElementList(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeElementList() start\n");
  
  if (_myRank  == _ioRank) {
    MPIIOFunctions::writeKeyValue<char>(fh, "\n!LIST_ELEM");
    MPIIOFunctions::writeKeyValue<char>(fh, "\n");
  }
  
  // get the local position in the file and broadcast it to all processes
  MPI_Offset offset;
  MPI_File_get_position(*fh, &offset);
  MPI_Bcast(&offset, 1, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
  // wOffset is initialized with current offset
  vector<MPI_Offset> wOffset(_nbWriters, offset); 
  
  // MeshElementType stores GLOBAL element type data
  const vector<MeshElementType>& me = MeshDataStack::getActive()->getTotalMeshElementTypes();
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = me.size();
  const CFuint nbLocalElements = getWriteData().getNbElements();
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);
  
  // store global info about the global ID ranges for sending
  // and the size of each send
  // each send will involve one writing process which wil collect all data 
  CFuint maxElemSendSize = 0;
  CFuint totalSize = 0;
  CFuint totalToSend = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbElementsInType = me[iType].elementCount;
    const CFuint nodesPlusStates  = me[iType].elementNodes + me[iType].elementStates;
    totalSize += nbElementsInType*nodesPlusStates;
    
    // fill in the writer list
    elementList.fill(nbElementsInType, nodesPlusStates, totalToSend);
    
    // update the maximum possible element-list size to send
    maxElemSendSize = max(maxElemSendSize, elementList.getMaxElemSize());
  }
  
  // sanity check
  cf_assert(totalToSend == totalSize);
  
  // start element list offset (current position)
  _offset[0].elems.first = offset;
  // end element list offset
  _offset[0].elems.second = _offset[0].elems.first + sizeof(CFuint)*totalSize;
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeElementList() => offsets = [" 
	<<  _offset[0].elems.first << ", " << _offset[0].elems.second << "]\n");
  
  Common::SafePtr< vector<CFuint> > globalElementIDs = 
    MeshDataStack::getActive()->getGlobalElementIDs();
  cf_assert(globalElementIDs->size() == nbLocalElements);
  
  CFLog(DEBUG_MAX,_myRank << " " << CFPrintContainer<vector<CFuint> >
		(" globalElementIDs  = ", &(*globalElementIDs)) << "\n");
  
  // ElementTypeData stores LOCAL element type data
  SafePtr< vector<ElementTypeData> > et = MeshDataStack::getActive()->getElementTypeData();
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  CFuint elemID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbLocalElementsInType = (*et)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem, ++elemID) {
      elementList.insertElemLocalID(elemID, (*globalElementIDs)[elemID], iType);
    }
  }
  elementList.endElemInsertion(_myRank);
  
  CFLog(VERBOSE,_myRank << " maxElemSendSize = " << maxElemSendSize << "\n");
  
  SafePtr<TopologicalRegionSet> elements =
    MeshDataStack::getActive()->getTrs("InnerCells");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  cf_assert(states.size() > 0);
  
  // buffer data to send
  vector<CFuint> sendElements(maxElemSendSize, 0);
  vector<CFuint> elementToPrint(maxElemSendSize, 0);
  
  vector<CFuint> countstates;
  if(_nbProc == 1 && states.size() == getWriteData().getNbElements()) { 
    countstates.resize(states.size(), (CFuint)0);
  }
  
  PE::Group& wg = PE::getGroup("Writers");
  CFint wRank = -1; 
  CFuint rangeID = 0;
  MPI_Offset dataSize = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbNodesInType  = me[iType].elementNodes;
    const CFuint nbStatesInType = me[iType].elementStates;
    const CFuint nodesPlusStates = nbNodesInType + nbStatesInType;
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
	  
	  CFuint isend = sendElemID*nodesPlusStates;
	  for (CFuint in = 0; in < nbNodesInType; ++in, ++isend) {
	    const CFuint localNodeID = elements->getNodeID(localElemID, in);
	    if (isend >= sendElements.size()) {
	      CFLogInfo(_myRank << " nbNodesInType = " << nbNodesInType
			<< " node isend = " << isend << " , size = " << sendElements.size() << "\n");
	      cf_assert(isend < sendElements.size());
	    }
	    
	    sendElements[isend] = nodes[localNodeID]->getGlobalID();
	  }
	  
	  for (CFuint in = 0; in < nbStatesInType; ++in, ++isend) {
	    const CFuint localStateID = elements->getStateID(localElemID, in);
	    if (isend >= sendElements.size()) {
	      CFLogInfo(_myRank << " state isend = " << isend << " , size = " << sendElements.size() << "\n");
	      cf_assert(isend < sendElements.size());
	    }
	    const CFuint sID = states[localStateID]->getGlobalID();
	    sendElements[isend] = sID;
	    // sanity check for cell-centered FV meshes here
	    if(_nbProc == 1 && states.size() == getWriteData().getNbElements()) {countstates[sID]++;}
	  }
	}
	
	cf_assert(eSize*nodesPlusStates <= elementList.getSendDataSize(rangeID));
      }
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		    (" sendElements  = ", &sendElements, nodesPlusStates) << "\n");
      
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

      
      /// another sanity check
      if(_nbProc == 1 && states.size() == getWriteData().getNbElements()) {
	static CFuint count = 0;
	count = 0;
	for (CFuint i = 0; i < sendSize; ++i) {
	  if (sendElements[i] != elementToPrint[i]) {
	    CFLog(INFO, "#"<< ++count << " => sendElements[" << i << "] = " << sendElements[i] << " != "  
		  << "elementToPrint[" << i << "] = " << elementToPrint[i] << "\n"); 
	    exit(1);
	  }
	}
      }
      ///
      
      
      // the offsets for all writers with send ID > current must be incremented  
      for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
	cf_assert(sendSize > 0);
	wOffset[iw] += sendSize*sizeof(CFuint);
	CFLog(DEBUG_MIN, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
      }
      
      //reset the all sendElement list to 0
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	sendElements[i] = 0;
      }
      
      // update the count element for the current element type
      countElem += elementList.getSendDataSize(rangeID)/nodesPlusStates;
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		    (" elementToPrint  = ", &elementToPrint, nodesPlusStates) << "\n");
      
    } // end sending loop
    
    if (_isWriterRank) { 
      CFLog(DEBUG_MIN, "ParCFmeshBinaryFileWriter::writeElementList() => P[" << _myRank 
	    << "] => offset = " << wOffset[wRank] << "\n");
      
      // each writer can now concurrently write all the collected data (related to one element type)
      cf_assert(wRank >= 0);
      MPIIOFunctions::writeAll("ParCFmeshBinaryFileWriter::writeElementList()", fh, wOffset[wRank], &elementToPrint[0], 
			       wSendSize, _maxBuffSize);
      
      // test  
      // cf_assert(elementToPrint.size() == wSendSize);
      // MPIIOFunctions::testBuff("ParCFmeshBinaryFileWriter::writeElementList()", wOffset[wRank], &elementToPrint[0], wSendSize);
    }
    
    // reset all the elements to print to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      elementToPrint[i] = 0;
    } 
    
    // reset the offset to the end of this type
    const CFuint nbElementsInType = me[iType].elementCount;
    dataSize += nbElementsInType*nodesPlusStates;
    wOffset.assign(wOffset.size(), _offset[0].elems.first + dataSize*sizeof(CFuint));
  } // end loop on iType
  
  if (_nbProc == 1 && states.size() == getWriteData().getNbElements()) {
    // here we check that, before writing them, state IDs are all pushed into the output buffer 
    for (CFuint i = 0; i < countstates.size(); ++i) {
      if (countstates[i] != 1) {
	CFLog(INFO, "ParCFmeshBinaryFileWriter::writeElementList() => ERROR: countstates[" << i << "] = " << countstates[i] << " != 1\n");
	exit(1);
      }
    }
  }
  
  CFLogInfo("Element written \n"); 
  
  if (_isWriterRank) {
    MPI_Barrier(wg.comm);
    MPI_File_seek(*fh, _offset[0].elems.second, MPI_SEEK_SET);
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeElementList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeTrsData(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeTrsData() start\n");
  
  vector<vector<CFuint> >&  trsInfo =
    MeshDataStack::getActive()->getTotalTRSInfo();
  
  const vector<std::string>& trsNames =
    MeshDataStack::getActive()->getTotalTRSNames();
  
  const CFuint nbTRSs = trsInfo.size();
  _offset[0].TRS.resize(nbTRSs);
  
  // CFLog(INFO, _myRank << " writes TRS starting from " << _offset[0].elems.second << "\n");
  
  if (_myRank == _ioRank) {
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_TRSs ", false, nbTRSs);
    CFLog(VERBOSE,"!NB_TRSs " << nbTRSs << "\n");
  }
  
  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    if (_myRank == _ioRank) {
      MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n!TRS_NAME ");
      MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, trsNames[iTRS]);
      CFLog(VERBOSE,"!TRS_NAME " << trsNames[iTRS] << "\n");
      
      const CFuint nbTRsInTRS = trsInfo[iTRS].size();
      MPIIOFunctions::MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!NB_TRs ", false, nbTRsInTRS);
      CFLog(VERBOSE,"!NB_TRs "   << nbTRsInTRS << "\n");
      
      MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n!NB_GEOM_ENTS ");
      MPI_File_write(*fh, &trsInfo[iTRS][0], (int)nbTRsInTRS, 
		     MPIStructDef::getMPIType(&trsInfo[iTRS][0]), &_status);
      CFLog(VERBOSE,CFPrintContainer<const vector<CFuint> >("!NB_GEOM_ENTS ", &trsInfo[iTRS], nbTRsInTRS));
      
      // AL: probably this geom_type info could go away ...
      MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n!GEOM_TYPE ");
      MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, CFGeoEnt::Convert::to_str(CFGeoEnt::FACE));
      CFLog(VERBOSE,"!GEOM_TYPE " << CFGeoEnt::FACE << "\n");
    }
    
    writeGeoList(iTRS, fh);
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeTrsData() end\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void ParCFmeshBinaryFileWriter::writeNodeList(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeNodeList() start\n");
  
  if (_isWriterRank) {
    MPI_File_seek(*fh, _offset[0].TRS.back().second, MPI_SEEK_SET);
  }
  
  if (_myRank == _ioRank) {
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n!LIST_NODE");
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n");
  }
  
  MPI_Offset offset;
  MPI_File_get_position(*fh, &offset);
  MPI_Bcast(&offset, 1, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
  // wOffset is initialized with current offset
  vector<MPI_Offset> wOffset(_nbWriters, offset); 
  
  PE::Group& wg = PE::getGroup("Writers");
  if (_isWriterRank) {
    MPI_Barrier(wg.comm);
  }
  
  const CFuint totNbNodes = MeshDataStack::getActive()->getTotalNodeCount();
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);

  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = 1; // just nodes
  const CFuint nbLocalElements = nodes.size();
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);
  
  // store global info about the global ID ranges for sending
  // and the size of each send
  const CFuint nbExtraNodalVars = getWriteData().getNbExtraNodalVars();
  getWriteData().prepareNodalExtraVars();
  
  const bool storePastNodes = getWriteData().storePastNodes();
  const vector<CFuint>& nodalExtraVarsStrides = *getWriteData().getExtraNodalVarStrides();
  const CFuint totalNbExtraNodalVars = std::accumulate(nodalExtraVarsStrides.begin(),
						       nodalExtraVarsStrides.end(), 0);
  CFuint nodesStride = dim + totalNbExtraNodalVars;
  if (storePastNodes) {nodesStride += dim;}
  
  // fill in the writer ist
  CFuint totalToSend = 0;
  elementList.fill(totNbNodes, nodesStride, totalToSend);
  
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = elementList.getMaxElemSize();
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  const CFuint nbLocalElementsInType = nodes.size();
  for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem) {
    elementList.insertElemLocalID(iElem, nodes[iElem]->getGlobalID(), 0);
  }
  elementList.endElemInsertion(_myRank);
  
  // start(current position) / end nodes list offset
  _offset[0].nodes.first  = offset;
  _offset[0].nodes.second = _offset[0].nodes.first + sizeof(CFreal)*totalToSend;
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeNodeList() => offsets = [" 
	<<  _offset[0].nodes.first << ", " << _offset[0].nodes.second << "]\n");
  
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
	const CFuint globalElemID = nodes[localElemID]->getGlobalID();
	const CFuint sendElemID = globalElemID - countElem;

	if (nodes[localElemID]->isParUpdatable()) {
	  CFuint isend = sendElemID*nodesStride;
	  for (CFuint in = 0; in < dim; ++in, ++isend) {
	    cf_assert(isend < sendElements.size());
	    sendElements[isend] = (*nodes[localElemID])[in]*refL;
	  }
	  
	  if (storePastNodes) {
	    const RealVector* pastNodesValues = getWriteData().getPastNode(localElemID);
	    cf_assert(pastNodesValues->size() == (*nodes[localElemID]).size());
	    for (CFuint in = 0; in < pastNodesValues->size(); ++in, ++isend) {
	      cf_assert(isend < sendElements.size());
	      sendElements[isend] = (*pastNodesValues)[in];
	    }
	  }
	  
	  if (nbExtraNodalVars > 0) {
	    const RealVector& extraNodalValues = getWriteData().getExtraNodalValues(localElemID);
	    cf_assert(extraNodalValues.size() == totalNbExtraNodalVars);
	    for (CFuint in = 0; in < totalNbExtraNodalVars; ++in, ++isend) {
	      cf_assert(isend < sendElements.size());
	      sendElements[isend] = extraNodalValues[in];
	    }
	  }
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
    MPI_Op_create((MPI_User_function *)cmpAndTakeMaxAbs2, 1, &myMpiOp);
    MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
	       MPIStructDef::getMPIType(&sendElements[0]), myMpiOp, wg.globalRanks[is], _comm);
    
    // MPI_Allreduce(&sendElements[0], &elementToPrint[0], maxElemSendSize, 
    //               MPIStructDef::getMPIType(&sendElements[0]), myMpiOp, _comm);
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" elementToPrint  = ", &elementToPrint, nodesStride) << "\n");
    
    // the offsets for all writers with send ID > current must be incremented  
    for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
      cf_assert(sendSize > 0);
      wOffset[iw] += sendSize*sizeof(CFreal);
      CFLog(DEBUG_MIN, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
    }
    
    // reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      sendElements[i] = 0;
    }
    
    // update the count element for the current element type
    countElem += elementList.getSendDataSize(rangeID)/nodesStride;
  }
  
  if (_isWriterRank) { 
    CFLog(DEBUG_MIN, "ParCFmeshBinaryFileWriter::writeNodeList() => P[" << _myRank << "] => offset = " << wOffset[wRank] << "\n");
    // each writer can now concurrently write all the collected data (related to one element type)
    cf_assert(wRank >= 0);
    // cout << _myRank << " writes " << wSendSize << "\n";
    // MPI_File_write_at_all(*fh, wOffset[wRank], &elementToPrint[0], (int)wSendSize, 
    // MPIStructDef::getMPIType(&elementToPrint[0]), &_status); 
    
    MPIIOFunctions::writeAll("ParCFmeshBinaryFileWriter::writeNodeList()", fh, wOffset[wRank], &elementToPrint[0], 
			     wSendSize, _maxBuffSize);
  }
  
  //reset the all sendElement list to 0
  for (CFuint i = 0; i < maxElemSendSize; ++i) {
    elementToPrint[i] = 0;
  }
  
  if (_isWriterRank) {
    MPI_Barrier(wg.comm);
    MPI_File_seek(*fh, _offset[0].nodes.second, MPI_SEEK_SET);
  }
  
  CFLogInfo("Nodes written \n");

  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeNodeList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeStateList(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeStateList() start\n");

  getWriteData().prepareStateExtraVars();

  if (_isWriterRank) {
    MPI_File_seek(*fh, _offset[0].nodes.second, MPI_SEEK_SET);
  }
  
  if (_myRank == _ioRank) {
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<CFuint>(fh, "\n!LIST_STATE ", false, getWriteData().isWithSolution());
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n");
  }
  
  PE::Group& wg = PE::getGroup("Writers");
  if (getWriteData().isWithSolution()){
    MPI_Offset offset;
    MPI_File_get_position(*fh, &offset);
    MPI_Bcast(&offset, 1, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
    // wOffset is initialized with current offset
    vector<MPI_Offset> wOffset(_nbWriters, offset); 
    
    if (_isWriterRank) {
      MPI_Barrier(wg.comm);
    }
    
    const CFuint totNbStates = MeshDataStack::getActive()->getTotalStateCount();
    const CFuint dim = PhysicalModelStack::getActive()->getNbEq();
    
    DataHandle < Framework::State*, Framework::GLOBAL > states =
      MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
    cf_assert(states.size() > 0);
    
    const CFuint nSend = _nbWriters;
    const CFuint nbElementTypes = 1; // just states
    const CFuint nbLocalElements = states.size();
    
    WriteListMap elementList;
    elementList.reserve(nbElementTypes, nSend, nbLocalElements);
    
    // store global info about the global ID ranges for sending
    // and the size of each send
    const CFuint nbExtraStateVars = getWriteData().getNbExtraStateVars();
    bool storePastStates = getWriteData().storePastStates();
    bool storeInterStates = getWriteData().storeInterStates();
    
    const vector<CFuint>& stateExtraVarsStrides = *getWriteData().getExtraStateVarStrides();
    const CFuint totalNbExtraStateVars = std::accumulate(stateExtraVarsStrides.begin(),
							 stateExtraVarsStrides.end(), 0);

    CFuint statesStride = dim + totalNbExtraStateVars;
    if (storePastStates)  {statesStride += dim;}
    if (storeInterStates) {statesStride += dim;}
    
    // fill in the writer ist
    CFuint totalToSend = 0;
    elementList.fill(totNbStates, statesStride, totalToSend);
    
    // update the maximum possible element-list size to send
    const CFuint maxElemSendSize = elementList.getMaxElemSize();
    
    // insert in the write list the local IDs of the elements
    // the range ID is automatically determined inside the WriteListMap
    const CFuint nbLocalElementsInType = states.size();
    for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem) {
      elementList.insertElemLocalID(iElem, states[iElem]->getGlobalID(), 0);
    }
    elementList.endElemInsertion(_myRank);
  
    // start(current position) / end nodes list offset
    _offset[0].states.first  = offset;
    _offset[0].states.second = _offset[0].states.first + sizeof(CFreal)*totalToSend;
    CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeStateList() => offsets = [" 
	  <<  _offset[0].states.first << ", " << _offset[0].states.second << "]\n");
    
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
	  const CFuint globalElemID = states[localElemID]->getGlobalID();
	  const CFuint sendElemID = globalElemID - countElem;
	  
	  if (states[localElemID]->isParUpdatable()) {
	    CFuint isend = sendElemID*statesStride;
	    for (CFuint in = 0; in < dim; ++in, ++isend) {
	      cf_assert(isend < sendElements.size());
	      sendElements[isend] = (*states[localElemID])[in];
	    }
	    
	    if(storePastStates){
	      const RealVector* pastStatesValues = getWriteData().getPastState(localElemID);
	      cf_assert(pastStatesValues->size() == (*states[localElemID]).size());
	      for (CFuint in = 0; in < pastStatesValues->size(); ++in, ++isend) {
		cf_assert(isend < sendElements.size());
		sendElements[isend] = (*pastStatesValues)[in];
	      }
	    }
	    
	    if(storeInterStates){
	      const RealVector* interStatesValues = getWriteData().getInterState(localElemID);
	      cf_assert(interStatesValues->size() == (*states[localElemID]).size());
	      for (CFuint in = 0; in < interStatesValues->size(); ++in, ++isend) {
		cf_assert(isend < sendElements.size());
		sendElements[isend] = (*interStatesValues)[in];
	      }
	    }
	    
	    if (nbExtraStateVars > 0) {
	      const RealVector& extraStateValues = getWriteData().getExtraStateValues(localElemID);
	      cf_assert(extraStateValues.size() == totalNbExtraStateVars);
	      for (CFuint in = 0; in < totalNbExtraStateVars; ++in, ++isend) {
		cf_assert(isend < sendElements.size());
		sendElements[isend] = extraStateValues[in];
	      }
	    }
	  }
	}
	
	cf_assert(eSize*statesStride <= elementList.getSendDataSize(rangeID));
      }

      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		    (" sendElements  = ", &sendElements, statesStride) << "\n");

      const CFuint sendSize = elementList.getSendDataSize(rangeID);
      cf_assert(sendSize <= sendElements.size());
      cf_assert(sendSize <= elementToPrint.size());
      
      // if the rank corresponds to a writing process, record the size to send for this writer
      if (_isWriterRank && wg.globalRanks[is] == _myRank) {
	wSendSize = sendSize; // this should be the total sendsize in the range
	wRank = is;
      }
      
      MPI_Op myMpiOp;
      MPI_Op_create((MPI_User_function *)cmpAndTakeMaxAbs2, 1, &myMpiOp);
      MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize, 
		 MPIStructDef::getMPIType(&sendElements[0]), 
		 myMpiOp, wg.globalRanks[is], _comm);
      // MPI_Allreduce(&sendElements[0], &elementToPrint[0], maxElemSendSize, 
      //               MPIStructDef::getMPIType(&sendElements[0]), myMpiOp, _comm);
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		    (" elementToPrint  = ", &elementToPrint, statesStride) << "\n");
      
      // the offsets for all writers with send ID > current must be incremented  
      for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
	cf_assert(sendSize > 0);
	wOffset[iw] += sendSize*sizeof(CFreal);
	CFLog(DEBUG_MIN, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
      }
      
      // reset the all sendElement list to 0
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	sendElements[i] = 0;
      }
      
      // update the count element for the current element type
      countElem += elementList.getSendDataSize(rangeID)/statesStride;
    }
    
    if (_isWriterRank) { 
      CFLog(DEBUG_MIN, "ParCFmeshBinaryFileWriter::writeStateList() => P[" << _myRank << "] => offset = " << wOffset[wRank] << "\n");
      // each writer can now concurrently write all the collected data (related to one element type)
      cf_assert(wRank >= 0);
     // CFLog(VERBOSE, _myRank << " writes " << wSendSize << " starting at " << wOffset[wRank] << " \n");
      // MPI_File_write_at_all(*fh, wOffset[wRank], &elementToPrint[0], (int)wSendSize,
      // MPIStructDef::getMPIType(&elementToPrint[0]), &_status); 
      
      MPIIOFunctions::writeAll("ParCFmeshBinaryFileWriter::writeStateList()", fh, wOffset[wRank], &elementToPrint[0], 
			       wSendSize, _maxBuffSize);
    }
    
    //reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      elementToPrint[i] = 0;
    }
  }
  
  if (_isWriterRank) {
    MPI_Barrier(wg.comm);
    MPI_File_seek(*fh, _offset[0].states.second, MPI_SEEK_SET);
  }
  
  CFLogInfo("States written \n");
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeStateList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeGeoList(CFuint iTRS, MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeGeoList() start\n");
  
  const vector<vector<CFuint> >&  trsInfo =
    MeshDataStack::getActive()->getTotalTRSInfo();
  
  SafePtr<vector<vector<vector<CFuint> > > > globalGeoIDS =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();
  
  const CFuint nbTRsInTRS = (*globalGeoIDS)[iTRS].size();
  CFLog(VERBOSE,"nbTRsInTRS = " << nbTRsInTRS << "\n");
  
  CFMat<CFuint> nbNodesStatesInTRGeoTmp(nbTRsInTRS, 2, static_cast<CFuint>(0));
  CFMat<CFuint> nbNodesStatesInTRGeo(nbTRsInTRS, 2, static_cast<CFuint>(0));
  
  const vector<std::string>& trsNames =
    MeshDataStack::getActive()->getTotalTRSNames();

  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(trsNames[iTRS]);
  cf_assert(trs->getNbTRs() == nbTRsInTRS);

  SafePtr<TopologicalRegionSet> elements =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const bool isFVMCC = (elements->getNbStatesInGeo(0) == 1);

  for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
    const CFuint nbGeosInLocalTR = (*trs)[iTR]->getLocalNbGeoEnts();
    CFLog(VERBOSE,"nbGeosInLocalTR " << iTR << " is " << nbGeosInLocalTR << "\n");

    // if there is at least one geometric entity, consider the first of them
    // and get its number of states and nodes
    if (nbGeosInLocalTR > 0) {
      CFuint maxNbNodesInTRGeo = 0;
      CFuint maxNbStatesInTRGeo = 0;
      const CFuint nbTRGeos = (*trs)[iTR]->getLocalNbGeoEnts();
      for (CFuint iGeo = 0; iGeo < nbTRGeos; ++iGeo) {
	maxNbNodesInTRGeo = max(maxNbNodesInTRGeo, (*trs)[iTR]->getNbNodesInGeo(iGeo));
	maxNbStatesInTRGeo = max(maxNbStatesInTRGeo, (*trs)[iTR]->getNbStatesInGeo(iGeo));
      }

      // in cell center FVM only the first state per TRS face must be considered
      // the second one is a ghost one that doesn't have to be written !!
      nbNodesStatesInTRGeoTmp(iTR,0) = maxNbNodesInTRGeo;
      nbNodesStatesInTRGeoTmp(iTR,1) = (isFVMCC) ? 1 : maxNbStatesInTRGeo;
    }
  }

  // fill in all the MAXIMUM numbers of nodes and states per all the TRs in this TRS
  // example: in a TR with quads and triangles nbNodesStatesInTRGeoTmp(iTR,0) = 4
  // example: in a TR with quads and triangles nbNodesStatesInTRGeoTmp(iTR,1) = 4 (FEM) or 1 (FVMCC)
  MPI_Allreduce(&nbNodesStatesInTRGeoTmp[0], &nbNodesStatesInTRGeo[0],
		nbNodesStatesInTRGeo.size(), MPIStructDef::getMPIType(&nbNodesStatesInTRGeoTmp[0]), 
		MPI_MAX, _comm);
  
  CFLog(VERBOSE,"nbNodesStatesInTRGeo = " << nbNodesStatesInTRGeo << "\n");
  
  if (_myRank == _ioRank) {
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n!LIST_GEOM_ENT ");
    // AL: this is a change to the old format: max number of nodes and states for B faces in TRS is written
    MPI_File_write(*fh, &nbNodesStatesInTRGeo[0], (int)nbNodesStatesInTRGeo.size(),
		   MPIStructDef::getMPIType(&nbNodesStatesInTRGeo[0]), &_status); 
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n");
  }
 
  MPI_Offset offset;
  MPI_File_get_position(*fh, &offset);
  MPI_Bcast(&offset, 1, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
  // wOffset is initialized with current offset
  vector<MPI_Offset> wOffset(_nbWriters, offset); 
  
  PE::Group& wg = PE::getGroup("Writers");
  if (_isWriterRank) {
    MPI_Barrier(wg.comm);
  }
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = nbTRsInTRS;
  const CFuint nbLocalElements = trs->getLocalNbGeoEnts();
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);

  // store global info about the global ID ranges for sending
  // and the size of each send
  CFuint maxElemSendSize = 0;
  CFuint totalSize = 0;
  CFuint totalToSend = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbElementsInType = trsInfo[iTRS][iType];
    const CFuint maxNbNodesInType =  nbNodesStatesInTRGeo(iType, 0);
    const CFuint maxNbStatesInType = nbNodesStatesInTRGeo(iType, 1);
    // the max number of nodes and states includes two additional integers:
    // one for the number of nodes, the second for the number of states in
    // the current TR geometric entity
    const CFuint maxNodesPlusStatesData = maxNbNodesInType + maxNbStatesInType + 2;
    totalSize += nbElementsInType*maxNodesPlusStatesData;
    
    // fill in the writer ist
    elementList.fill(nbElementsInType, maxNodesPlusStatesData, totalToSend);
    
    // update the maximum possible element-list size to send
    maxElemSendSize = max(maxElemSendSize, elementList.getMaxElemSize());
  }
  
  // sanity check
  cf_assert(totalToSend == totalSize);
  
  // start TRS element list offset (current position)
  _offset[0].TRS[iTRS].first  = offset;
  // end TRS element list offset
  _offset[0].TRS[iTRS].second = _offset[0].TRS[iTRS].first + sizeof(CFint)*totalToSend;
  
  // cout << _myRank << " TRS[" << iTRS << "] offsets = [" << _offset[0].TRS[iTRS].first  << ", " << _offset[0].TRS[iTRS].second << "]\n";
  
  if (_isWriterRank) {
    MPI_File_seek(*fh, _offset[0].TRS[iTRS].first, MPI_SEEK_SET);
  }
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbLocalElementsInType = (*trs)[iType]->getLocalNbGeoEnts();
    for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem) {
      elementList.insertElemLocalID(iElem, (*globalGeoIDS)[iTRS][iType][iElem], iType);
    }
  }
  elementList.endElemInsertion(_myRank);
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  cf_assert(states.size() > 0);
  
  // buffer data to send
  vector<CFint> sendElements(maxElemSendSize, -1);
  vector<CFint> elementToPrint(maxElemSendSize, -1);
  
  CFint wRank = -1; 
  CFuint rangeID = 0;
  MPI_Offset dataSize = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint maxNbNodesInType  = nbNodesStatesInTRGeo(iType, 0);
    const CFuint maxNbStatesInType = nbNodesStatesInTRGeo(iType, 1);
    const CFuint maxNodesPlusStates = maxNbNodesInType + maxNbStatesInType + 2;
    
    CFuint wSendSize = 0;
    CFuint countElem = 0;
    for (CFuint is = 0; is < nSend; ++is, ++rangeID) {
      bool isRangeFound = false;
      WriteListMap::List elist = elementList.find(rangeID, isRangeFound);
      
      if (isRangeFound) {
	CFuint eSize = 0;
	for (WriteListMap::ListIterator it = elist.first; it != elist.second; ++it, ++eSize) {
	  const CFuint localElemID = it->second;
	  const CFuint globalElemID = (*globalGeoIDS)[iTRS][iType][localElemID];
	  const CFuint sendElemID = globalElemID - countElem;
	  
	  CFLogDebugMax("localElemID = " << localElemID <<
			", globalElemID = " << globalElemID <<
			", sendElemID = " << sendElemID << "\n");
	  
	  CFuint isend = sendElemID*maxNodesPlusStates;
	  // number of nodes in the current TR geo entity
	  const CFuint nbNodesInTRGeo  = (*trs)[iType]->getNbNodesInGeo(localElemID);
	  sendElements[isend++] = nbNodesInTRGeo;
	  
	  // number of states in the current TR geo entity
	  const CFuint nbStatesInTRGeo = (isFVMCC) ? 1 : (*trs)[iType]->getNbStatesInGeo(localElemID);
	  sendElements[isend++] = nbStatesInTRGeo;

	  // TR geo nodes data
	  for (CFuint in = 0; in < maxNbNodesInType; ++in, ++isend) {
	    // the local node ID is set to -1 if the maximum number of nodes exceeds the actual value
	    const CFint localNodeID = (in < nbNodesInTRGeo) ?
	      static_cast<CFint>((*trs)[iType]->getNodeID(localElemID, in)) : -1;
	    if (isend >= sendElements.size()) {
	      CFLog(VERBOSE,_myRank << " maxNbNodesInType = " << maxNbNodesInType
			    << " node isend = " << isend << " , size = " << sendElements.size() << "\n");
	      cf_assert(isend < sendElements.size());
	    }
	    // set the global ID to -1 if the local ID is -1
	    sendElements[isend] = (localNodeID != -1) ?
	      static_cast<CFint>(nodes[localNodeID]->getGlobalID()) : -1;
	  }
	  
	  // TR geo states data
	  for (CFuint in = 0; in < maxNbStatesInType; ++in, ++isend) {
	    // the local state ID is set to -1 if the maximum number of states exceeds the actual value
	    const CFint localStateID = (in < nbStatesInTRGeo) ?
	      static_cast<CFint>((*trs)[iType]->getStateID(localElemID, in)) : -1;
	    if (isend >= sendElements.size()) {
	      CFLog(VERBOSE,_myRank << " state isend = " << isend << " , size = " << sendElements.size() << "\n");
	      cf_assert(isend < sendElements.size());
	    }
	    // set the global ID to -1 if the local ID is -1
	    sendElements[isend] = (localStateID != -1) ?
	      static_cast<CFint>(states[localStateID]->getGlobalID()) : -1;
	  }
	}
	
	cf_assert(eSize*maxNodesPlusStates <= elementList.getSendDataSize(rangeID));
      }
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFint> >
		    (" sendElements  = ", &sendElements, maxNodesPlusStates) << "\n");
      
      const CFuint sendSize = elementList.getSendDataSize(rangeID);
      cf_assert(sendSize <= sendElements.size());
      cf_assert(sendSize <= elementToPrint.size());
      
      // if the rank corresponds to a writing process, record the size to send for this writer
      if (_isWriterRank && wg.globalRanks[is] == _myRank) {
	wSendSize = sendSize; // this should be the total sendsize in the range
	wRank = is;
      }
      
      // for each send, accumulate data to the corresponding writing process
      MPI_Reduce(&sendElements[0], &elementToPrint[0], (int)sendSize,
		 MPIStructDef::getMPIType(&sendElements[0]), MPI_MAX, wg.globalRanks[is], _comm);
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFint> >
		    (" elementToPrint  = ", &elementToPrint, maxNodesPlusStates) << "\n");
      
      // the offsets for all writers with send ID > current must be incremented  
      for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
	cf_assert(sendSize > 0);
	wOffset[iw] += sendSize*sizeof(CFint);
	CFLog(DEBUG_MIN, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
      }
      
      //reset the all sendElement list to -1
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	sendElements[i] = -1;
      }
      
      // update the count element for the current element type
      countElem += elementList.getSendDataSize(rangeID)/maxNodesPlusStates;
    } // end sending loop
    
    if (_isWriterRank) {
      CFLog(DEBUG_MIN, "ParCFmeshBinaryFileWriter::writeGeoList() => P[" << _myRank 
	    << "] => offset = " << wOffset[wRank] << "\n");
      cf_assert(wRank >= 0);
      // CFLog(VERBOSE, "P[" << _myRank << "] writes " << wSendSize << " starting from " << wOffset[wRank] << "\n");
      // MPI_File_write_at_all(*fh, wOffset[wRank], &elementToPrint[0], (int)wSendSize, 
      //MPIStructDef::getMPIType(&elementToPrint[0]), &_status); 
      
      MPIIOFunctions::writeAll("ParCFmeshBinaryFileWriter::writeGeoList()", fh, wOffset[wRank], &elementToPrint[0], 
			       wSendSize, _maxBuffSize);
      
      // note here we are writing some components that could be = -1
      // this will have to be taken into account in the parallel reader 
    }
    
    // reset all the elements to print to -1
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      elementToPrint[i] = -1;
    }
    
    // reset the offset to the end of this type
    const CFuint nbElementsInType = trsInfo[iTRS][iType];
    dataSize += nbElementsInType*maxNodesPlusStates;
    wOffset.assign(wOffset.size(), _offset[0].TRS[iTRS].first + dataSize*sizeof(CFint));
    CFLog(DEBUG_MIN, CFPrintContainer<vector<MPI_Offset> >(" TR wOffset = ", &wOffset));
  }
  
  if (_isWriterRank) {
    MPI_Barrier(wg.comm);
    MPI_File_seek(*fh, _offset[0].TRS[iTRS].second, MPI_SEEK_SET);
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeGeoList() end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileWriter::writeEndFile(MPI_File* fh)
{
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeEndFile() start\n");
  
  if (_myRank == _ioRank) {
    MPI_Offset offset;
    MPI_File_get_position(*fh, &offset);
    CFLog(VERBOSE, _myRank << " END FILE offset is " << offset << "\n");
    MPIIOFunctions::MPIIOFunctions::writeKeyValue<char>(fh, "\n!END");
  }
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeEndFile() end\n");
}

//////////////////////////////////////////////////////////////////////////////
 
    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
