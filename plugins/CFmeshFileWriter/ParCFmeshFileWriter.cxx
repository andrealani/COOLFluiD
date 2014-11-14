// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>
#include <numeric>

#include "ParCFmeshFileWriter.hh"
#include "Framework/ElementTypeData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/ElementTypeData.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Common/CFMultiMap.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/MPI/MPIStructDef.hh"

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

void cmpAndTakeMaxAbs(CFreal* invec, CFreal* inoutvec, int* len,
		      MPI_Datatype* datatype)
{
  cf_assert(len != CFNULL);
  int size = *len;
  for (int i = 0; i < size; ++i) {
    inoutvec[i] = (fabs(invec[i]) > 0.) ? invec[i] : inoutvec[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

ParCFmeshFileWriter::ParCFmeshFileWriter() :
  ConfigObject("ParCFmeshFileWriter"),
  _comm(0),
  _myRank(0),
  _nbProc(0),
  _ioRank(0),
  _isNewFile(0),
  _writeData(),
  _fileList(),
  _mapFileToStartNodeList()
{
  addConfigOptionsTo(this);
  
  _comm   = PE::GetPE().GetCommunicator();
  _myRank = PE::GetPE().GetRank();
  _nbProc = PE::GetPE().GetProcessorCount();
  cf_assert(_nbProc > 0);
}

//////////////////////////////////////////////////////////////////////////////

ParCFmeshFileWriter::~ParCFmeshFileWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeToFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  ofstream* file = CFNULL;

  // reset to 0 the new file flag
  _isNewFile = 0;

   // AL: the processor with less elements is chosen as the master IO node
   // this can avoid memory "explosion" on master node in cases for which 
   // the partitioning is not well balanced
   CFuint nbLocalElements = getWriteData().getNbElements();
   CFuint minNumberElements = 0;
   MPI_Allreduce(&nbLocalElements, &minNumberElements, 1, MPIStructDef::getMPIType(&nbLocalElements), MPI_MIN, _comm);
   CFuint rank = (minNumberElements == nbLocalElements)  ? _myRank : 0;
   // IO rank is maximum rank whose corresponding process has minimum number of elements
   MPI_Allreduce(&rank, &_ioRank, 1, MPIStructDef::getMPIType(&rank), MPI_MAX, _comm);    
   CFout << "ParCFmeshFileWriter::writeToFile() => IO rank is " << _ioRank << "\n";
 
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  if (_myRank == _ioRank) {
    // if the file has already been processed once, open in I/O mode
    if (_fileList.count(filepath) > 0) {
      file = &fhandle->open(filepath, ios_base::in | ios_base::out);
    }
    else {
      file = &fhandle->open(filepath);

      // if the file is a new one add it to the file list
      _fileList.insert(filepath);
      _isNewFile = 1;
    }
  }

  // communicate to all the processors about the status of the current file
  MPI_Bcast(&_isNewFile, 1, MPIStructDef::getMPIType(&_isNewFile), _ioRank, _comm);

  writeToFileStream(filepath, file);

  if (_myRank == _ioRank) {
    fhandle->close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeToFileStream
(const boost::filesystem::path& filepath,
 std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeToFileStream() start\n");

  if (_isNewFile == 1)
  {
    writeVersionStamp(fout);
    // global counts
    writeGlobalCounts(fout);

    // extra vars info
    writeExtraVarsInfo(fout);

    // elements info
    writeElements(fout);

    // TRS data
    writeTrsData(fout);

    if (_myRank == _ioRank) {
      _mapFileToStartNodeList[filepath] = fout->tellp();
    }
  }

  // what follows works but can leave some lines at the end of the file
  // due to previous writing ...
  else {
    cf_assert(_isNewFile == 0);

    long position = 0;
    if (_myRank == _ioRank) {
      position = _mapFileToStartNodeList.find(filepath)->second;
      fout->seekp(position, ios::beg);
    }

    // communicate to all the processors about the position in the current file
    MPI_Bcast(&position, 1, MPI_LONG, _ioRank, _comm);
  }

  
  // extra vars that are not state or node related
  writeExtraVars(fout);
  
  // write the node list
  writeNodeList(fout);

  // write the state list
  writeStateList(fout);

  if (_myRank == _ioRank) {
    *fout << "!END" << "\n";
  }

  CFLogDebugMin( "ParCFmeshFileWriter::writeFile() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeVersionStamp(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeDimension() called" << "\n");

  if (_myRank  == _ioRank)
  {
    *fout << "!COOLFLUID_VERSION "    << CFEnv::getInstance().getCFVersion() << "\n";
    *fout << "!COOLFLUID_SVNVERSION " << CFEnv::getInstance().getSvnVersion() << "\n";
    *fout << "!CFMESH_FORMAT_VERSION 1.3\n";
  }

  CFLogDebugMin( "ParCFmeshFileWriter::writeDimension() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeGlobalCounts(std::ofstream *const fout)
{
 CFLogDebugMin( "ParCFmeshFileWriter::writeDimension() called" << "\n");

 if (_myRank  == _ioRank) {
   *fout << "!NB_DIM " << PhysicalModelStack::getActive()->getDim() << "\n";

   *fout << "!NB_EQ " << PhysicalModelStack::getActive()->getNbEq() << "\n";

   *fout << "!NB_NODES " << MeshDataStack::getActive()->getTotalNodeCount() << " "
	 << getWriteData().getNbNonUpdatableNodes() << "\n";

   *fout << "!NB_STATES " << MeshDataStack::getActive()->getTotalStateCount() << " "
	 << getWriteData().getNbNonUpdatableStates() << "\n";

   const std::vector<CFuint>& tElem =  MeshDataStack::getActive()->getTotalElements();
   *fout << "!NB_ELEM " << std::accumulate(tElem.begin(), tElem.end(), 0) << "\n";
 }

  CFLogDebugMin( "ParCFmeshFileWriter::writeDimension() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeExtraVarsInfo(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeExtraVarsInfo() called" << "\n");

 if (_myRank  == _ioRank) {

   if(getWriteData().storePastStates()) {
     *fout << "!STORE_PASTSTATES " << getWriteData().storePastStates() << "\n";
   }

   if(getWriteData().storePastNodes()) {
     *fout << "!STORE_PASTNODES " << getWriteData().storePastNodes() << "\n";
   }

  if(getWriteData().storeInterStates()) {
     *fout << "!STORE_INTERSTATES " << getWriteData().storeInterStates() << "\n";
   }

   if(getWriteData().storeInterNodes()) {
     *fout << "!STORE_INTERNODES " << getWriteData().storeInterNodes() << "\n";
   }

   const CFuint nbExtraNodalVars = getWriteData().getNbExtraNodalVars();
   if(nbExtraNodalVars > 0) {
     *fout << "!NB_EXTRA_NVARS " << nbExtraNodalVars << "\n";

     *fout << "!EXTRA_NVARS_NAMES ";
     for(CFuint iVar = 0; iVar < nbExtraNodalVars; iVar++) {
       *fout << (*(getWriteData().getExtraNodalVarNames()))[iVar] << " ";
     }
     *fout << "\n";

     *fout << "!EXTRA_NVARS_STRIDES ";
     for(CFuint iVar = 0; iVar < nbExtraNodalVars; iVar++) {
       *fout << (*(getWriteData().getExtraNodalVarStrides()))[iVar] << " ";
     }
     *fout << "\n";
   }

   const CFuint nbExtraStateVars = getWriteData().getNbExtraStateVars();
   if (nbExtraStateVars > 0) {
     *fout << "!NB_EXTRA_SVARS " << nbExtraStateVars << "\n";

     *fout << "!EXTRA_SVARS_NAMES ";
     for(CFuint iVar = 0; iVar < nbExtraStateVars; iVar++) {
       *fout << (*(getWriteData().getExtraStateVarNames()))[iVar] << " ";
     }
     *fout << "\n";

     *fout << "!EXTRA_SVARS_STRIDES ";
     for(CFuint iVar = 0; iVar < nbExtraStateVars; iVar++) {
       *fout << (*(getWriteData().getExtraStateVarStrides()))[iVar] << " ";
     }
     *fout << "\n";
   }
   
   const CFuint nbExtraVars = getWriteData().getNbExtraVars();
   if (nbExtraVars > 0) {
      *fout << "!NB_EXTRA_VARS " << nbExtraVars << "\n";
   
      *fout << "!EXTRA_VARS_NAMES ";
      for(CFuint iVar = 0; iVar < nbExtraVars; iVar++) {
        *fout << (*(getWriteData().getExtraVarNames()))[iVar] << " ";
      }
      *fout << "\n";
      
      *fout << "!EXTRA_VARS_STRIDES ";
      for(CFuint iVar = 0; iVar < nbExtraVars; iVar++) {
        *fout << (*(getWriteData().getExtraVarStrides()))[iVar] << " ";
      }
      *fout << "\n";
    }
   
 }

  CFLogDebugMin( "ParCFmeshFileWriter::writeExtraVarsInfo() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////
      

void ParCFmeshFileWriter::writeExtraVars(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeNodeList() called" << "\n");
  
  if (_myRank == _ioRank) {
    
    getWriteData().prepareExtraVars();
    const vector<CFuint>& extraVarsStrides = *getWriteData().getExtraVarStrides();
    const CFuint totalNbExtraVars = std::accumulate(extraVarsStrides.begin(),extraVarsStrides.end(), 0);
    const RealVector& extraValues = getWriteData().getExtraValues();
    cf_assert(extraValues.size() == totalNbExtraVars);
    
    cf_assert(fout != CFNULL);
    *fout << "!EXTRA_VARS ";
    for(CFuint iVar = 0; iVar < totalNbExtraVars; iVar++) {
      *fout << extraValues[iVar] << " ";
    }
    *fout << "\n";
  }
}  

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeElements(std::ofstream *const fout)
{
  CFLogDebugMin("ParCFmeshFileWriter::writeElements() start\n");

  if (_myRank  == _ioRank) {
    *fout << "!NB_ELEM_TYPES " << MeshDataStack::getActive()->getTotalElements().size() << "\n";
    *fout << "!GEOM_POLYORDER " << static_cast<CFuint> (getWriteData().getGeometricPolyOrder()) << "\n";
    *fout << "!SOL_POLYORDER "  << static_cast<CFuint> (getWriteData().getSolutionPolyOrder()) << "\n";
    *fout << "!ELEM_TYPES ";

    const vector<MeshElementType>& me =
      MeshDataStack::getActive()->getTotalMeshElementTypes();

    const CFuint nbElementTypes = me.size();
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      *fout << me[i].elementName << " ";
    }
    *fout << "\n";

    *fout << "!NB_ELEM_PER_TYPE";
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      *fout << " " << me[i].elementCount;
    }
    *fout << "\n";

    *fout << "!NB_NODES_PER_TYPE";
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      *fout << " " << me[i].elementNodes;
    }
    *fout << "\n";

    *fout << "!NB_STATES_PER_TYPE";
    for (CFuint i = 0; i < nbElementTypes; ++i) {
      *fout << " " << me[i].elementStates;
    }
    *fout << "\n";
  }

  writeElementList(fout);

  CFLogDebugMin( "ParCFmeshFileWriter::writeElements() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeElementList(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeElementList() called" << "\n");

  if (_myRank  == _ioRank) {
    *fout << "!LIST_ELEM " << "\n";
  }

  // MeshElementType stores GLOBAL element type data
  const vector<MeshElementType>& me =
    MeshDataStack::getActive()->getTotalMeshElementTypes();

  const CFuint nSend = _nbProc;
  const CFuint nbElementTypes = me.size();
  const CFuint nbLocalElements = getWriteData().getNbElements();

  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);

  // store global info about the global ID ranges for sending
  // and the size of each send
  CFuint maxElemSendSize = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbElementsInType = me[iType].elementCount;
    const CFuint nbNodesInType =  me[iType].elementNodes;
    const CFuint nbStatesInType =  me[iType].elementStates;
    const CFuint nodesPlusStates = nbNodesInType + nbStatesInType;
    const CFuint maxElemSize = (nbElementsInType/nSend + nbElementsInType%nSend)*nodesPlusStates;
    const CFuint minElemSize = (nbElementsInType/nSend)*nodesPlusStates;

    // update the maximum possible element-list size to send
    maxElemSendSize = max(maxElemSendSize, maxElemSize);

    CFuint minElemID = 0;
    for (CFuint is = 0; is < nSend; ++is) {
      CFuint currSendSize = 0;
      if (is == 0) {
	// first send per iType has maximum size
	currSendSize = maxElemSize;

      }
      else {
	// all the other sends (is != 0) per iType have minimum size
	currSendSize = minElemSize;
      }
      // add the global size of the element data (each one of size=nodesPlusStates)
      // to communicate during this send
      elementList.addSendDataSize(currSendSize);
      const CFuint maxElemID = minElemID + currSendSize/nodesPlusStates;

      // insert the range
      elementList.addRange(minElemID, maxElemID);

      // set the new minimum as the old maximum
      minElemID = maxElemID;
    }
  }

  Common::SafePtr< vector<CFuint> > globalElementIDs =
    MeshDataStack::getActive()->getGlobalElementIDs();
  cf_assert(globalElementIDs->size() == nbLocalElements);

  CFLogDebugMin(_myRank << " " << CFPrintContainer<vector<CFuint> >
		(" globalElementIDs  = ", &(*globalElementIDs)) << "\n");

  // ElementTypeData stores LOCAL element type data
  SafePtr< vector<ElementTypeData> > et =
    MeshDataStack::getActive()->getElementTypeData();

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

  CFLogDebugMin(_myRank << " maxElemSendSize = " << maxElemSendSize << "\n");

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

  CFuint rangeID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbNodesInType  = me[iType].elementNodes;
    const CFuint nbStatesInType = me[iType].elementStates;
    const CFuint nodesPlusStates = nbNodesInType + nbStatesInType;

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
	    sendElements[isend] = states[localStateID]->getGlobalID();
	  }
	}

	cf_assert(eSize*nodesPlusStates <= elementList.getSendDataSize(rangeID));
      }

    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
	    (" sendElements  = ", &sendElements, nodesPlusStates) << "\n");

      const CFuint sendSize = elementList.getSendDataSize(rangeID);
      cf_assert(sendSize <= sendElements.size());
      cf_assert(sendSize <= elementToPrint.size());

//      MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
//		 MPIStructDef::getMPIType(&sendElements[0]), MPI_MAX, _ioRank, _comm); // allocates too much memory on master node with OPENMPI
      
      MPI_Allreduce(&sendElements[0], &elementToPrint[0], maxElemSendSize, MPIStructDef::getMPIType(&sendElements[0]), MPI_MAX, _comm);
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		    (" elementToPrint  = ", &elementToPrint, nodesPlusStates) << "\n");

      if (_myRank == _ioRank) {
	for (CFuint i = 0; i < sendSize; ++i) {
	  if ((i+1)%nodesPlusStates > 0) {
	    *fout << elementToPrint[i] << " ";
	  }
	  else {
	    *fout << elementToPrint[i] << endl;
	  }
	}
      }

      //reset the all sendElement list to 0
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	sendElements[i] = 0;
	elementToPrint[i] = 0;
      }

      // update the count element for the current element type
      countElem += elementList.getSendDataSize(rangeID)/nodesPlusStates;
    }
  }
  
  CFLogInfo("Element written \n"); 

  CFLogDebugMin( "ParCFmeshFileWriter::writeElementList() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeTrsData(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeTrsData() called" << "\n");

  const vector<vector<CFuint> >&  trsInfo =
    MeshDataStack::getActive()->getTotalTRSInfo();

  const vector<std::string>& trsNames =
    MeshDataStack::getActive()->getTotalTRSNames();

  const CFuint nbTRSs = trsInfo.size();
  if (_myRank == _ioRank) {
    *fout << "!NB_TRSs " << nbTRSs << "\n";
    CFLogDebugMin("!NB_TRSs " << nbTRSs << "\n");
  }

  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    if (_myRank == _ioRank) {
      *fout << "!TRS_NAME " << trsNames[iTRS] << "\n";
      CFLogDebugMin("!TRS_NAME " << trsNames[iTRS] << "\n");

      const CFuint nbTRsInTRS = trsInfo[iTRS].size();
      *fout << "!NB_TRs "   << nbTRsInTRS << "\n";
      CFLogDebugMin("!NB_TRs "   << nbTRsInTRS << "\n");

      *fout << "!NB_GEOM_ENTS";
      for(CFuint tr = 0; tr < nbTRsInTRS; ++tr) {
	*fout << " " << trsInfo[iTRS][tr];
      }
      *fout << "\n";

      CFLogDebugMin(CFPrintContainer<const vector<CFuint> >("!NB_GEOM_ENTS ", &trsInfo[iTRS], nbTRsInTRS));

      // AL: probably this geom_type info could go away ...
      *fout << "!GEOM_TYPE " << CFGeoEnt::FACE << "\n";

      CFLogDebugMin("!GEOM_TYPE " << CFGeoEnt::FACE << "\n");
    }

    writeGeoList(iTRS, fout);
  }

  CFLogDebugMin( "ParCFmeshFileWriter::writeTrsData() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeNodeList(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeNodeList() called" << "\n");

  if (_myRank == _ioRank) {
    cf_assert(fout != CFNULL);
    *fout << "!LIST_NODE " << "\n";
  }

  const CFuint totNbNodes = MeshDataStack::getActive()->getTotalNodeCount();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal refL = PhysicalModelStack::getActive()->
    getImplementor()->getRefLength();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);

  const CFuint nSend = _nbProc;
  const CFuint nbElementTypes = 1; // just nodes
  const CFuint nbLocalElements = nodes.size();

  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);

  // store global info about the global ID ranges for sending
  // and the size of each send
  const CFuint nbExtraNodalVars = getWriteData().getNbExtraNodalVars();
  getWriteData().prepareNodalExtraVars();

  bool storePastNodes = getWriteData().storePastNodes();

  const vector<CFuint>& nodalExtraVarsStrides = *getWriteData().getExtraNodalVarStrides();
  const CFuint totalNbExtraNodalVars = std::accumulate(nodalExtraVarsStrides.begin(),
						  nodalExtraVarsStrides.end(), 0);
  CFuint nodesStride = dim + totalNbExtraNodalVars;

  if (storePastNodes) {
    nodesStride += dim;
  }

  const CFuint maxElemSize = (totNbNodes/nSend + totNbNodes%nSend)*nodesStride;
  const CFuint minElemSize = (totNbNodes/nSend)*nodesStride;
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = maxElemSize;

  CFuint minElemID = 0;
  for (CFuint is = 0; is < nSend; ++is) {
    CFuint currSendSize = 0;
    if (is == 0) {
      // first send per iType has maximum size
      currSendSize = maxElemSize;
    }
    else {
      // all the other sends (is != 0) per iType have minimum size
      currSendSize = minElemSize;
    }
    // add the global size of the element data (each one of size=nodesPlusStates)
    // to communicate during this send
    elementList.addSendDataSize(currSendSize);
    const CFuint maxElemID = minElemID + currSendSize/nodesStride;

    // insert the range
    elementList.addRange(minElemID, maxElemID);

    // set the new minimum as the old maximum
    minElemID = maxElemID;
  }

  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  const CFuint nbLocalElementsInType = nodes.size();
  for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem) {
    elementList.insertElemLocalID(iElem, nodes[iElem]->getGlobalID(), 0);
  }
  elementList.endElemInsertion(_myRank);

  // buffer data to send
  vector<CFreal> sendElements(maxElemSendSize, 0);
  vector<CFreal> elementToPrint(maxElemSendSize, 0);

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

	CFuint isend = sendElemID*nodesStride;
	for (CFuint in = 0; in < dim; ++in, ++isend) {
	  cf_assert(isend < sendElements.size());
	  sendElements[isend] = (*nodes[localElemID])[in];
	}

        if(storePastNodes){
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

      cf_assert(eSize*nodesStride <= elementList.getSendDataSize(rangeID));
    }

    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" sendElements  = ", &sendElements, nodesStride) << "\n");

    const CFuint sendSize = elementList.getSendDataSize(rangeID);
    cf_assert(sendSize <= sendElements.size());
    cf_assert(sendSize <= elementToPrint.size());

    MPI_Op myMpiOp;
    MPI_Op_create((MPI_User_function *)cmpAndTakeMaxAbs, 1, &myMpiOp);
    
    MPI_Datatype MPI_CFREAL = MPIStructDef::getMPIType(&elementToPrint[0]);
    
//    MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
//	       MPI_CFREAL, myMpiOp, _ioRank, _comm);
    
    MPI_Allreduce(&sendElements[0], &elementToPrint[0], maxElemSendSize, MPI_CFREAL, myMpiOp, _comm);
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" elementToPrint  = ", &elementToPrint, nodesStride) << "\n");

    if (_myRank == _ioRank) {
      CFuint countN = 0;
      for (CFuint i = 0; i < sendSize; ++i) {
	const CFreal coeff = (countN < dim) ? refL : 1.;
	if ((i+1)%nodesStride > 0) {
          fout->precision(14);
          fout->setf(ios::scientific,ios::floatfield);
	  *fout << elementToPrint[i]*coeff << " ";
	  countN++;
	}
	else {
          fout->precision(14);
          fout->setf(ios::scientific,ios::floatfield);
	  *fout << elementToPrint[i]*coeff << endl;
	  countN = 0; // reset countN to 0
	}
      }
    }

    //reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      sendElements[i] = 0;
      elementToPrint[i] = 0;
    }

    // update the count element for the current element type
    countElem += elementList.getSendDataSize(rangeID)/nodesStride;
  }

  CFLogInfo("Nodes written \n");

  CFLogDebugMin( "ParCFmeshFileWriter::writeNodeList() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeStateList(std::ofstream *const fout)
{
  CFLogDebugMin( "ParCFmeshFileWriter::writeStateList() called" << "\n");

  getWriteData().prepareStateExtraVars();

  if (_myRank == _ioRank) {
    *fout << "!LIST_STATE " << getWriteData().isWithSolution() << "\n";
  }

  if (getWriteData().isWithSolution()){
    const CFuint totNbStates = MeshDataStack::getActive()->getTotalStateCount();
    const CFuint dim = PhysicalModelStack::getActive()->getNbEq();

    DataHandle < Framework::State*, Framework::GLOBAL > states =
      MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
    cf_assert(states.size() > 0);

    const CFuint nSend = _nbProc;
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

    if(storePastStates) {
      statesStride += dim;
    }
    if(storeInterStates) {
      statesStride += dim;
    }

    const CFuint maxElemSize = (totNbStates/nSend + totNbStates%nSend)*statesStride;
    const CFuint minElemSize = (totNbStates/nSend)*statesStride;
    // update the maximum possible element-list size to send
    const CFuint maxElemSendSize = maxElemSize;

    CFuint minElemID = 0;
    for (CFuint is = 0; is < nSend; ++is) {
      CFuint currSendSize = 0;
      if (is == 0) {
	// first send per iType has maximum size
	currSendSize = maxElemSize;
      }
      else {
	// all the other sends (is != 0) per iType have minimum size
	currSendSize = minElemSize;
      }
      // add the global size of the element data (each one of size=statesStride)
      // to communicate during this send
      elementList.addSendDataSize(currSendSize);
      const CFuint maxElemID = minElemID + currSendSize/statesStride;

      // insert the range
      elementList.addRange(minElemID, maxElemID);

      // set the new minimum as the old maximum
      minElemID = maxElemID;
    }

    // insert in the write list the local IDs of the elements
    // the range ID is automatically determined inside the WriteListMap
    const CFuint nbLocalElementsInType = states.size();
    for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem) {
      elementList.insertElemLocalID(iElem, states[iElem]->getGlobalID(), 0);
    }
    elementList.endElemInsertion(_myRank);

    // buffer data to send
    vector<CFreal> sendElements(maxElemSendSize, 0);
    vector<CFreal> elementToPrint(maxElemSendSize, 0);

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

	cf_assert(eSize*statesStride <= elementList.getSendDataSize(rangeID));
      }

      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		    (" sendElements  = ", &sendElements, statesStride) << "\n");

      const CFuint sendSize = elementList.getSendDataSize(rangeID);
      cf_assert(sendSize <= sendElements.size());
      cf_assert(sendSize <= elementToPrint.size());

      MPI_Op myMpiOp;
      MPI_Op_create((MPI_User_function *)cmpAndTakeMaxAbs, 1, &myMpiOp);

      MPI_Datatype MPI_CFREAL = MPIStructDef::getMPIType(&elementToPrint[0]);

      //  MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
      // 		 MPI_CFREAL, myMpiOp, _ioRank, _comm);
      
      MPI_Allreduce(&sendElements[0], &elementToPrint[0], maxElemSendSize, MPI_CFREAL, myMpiOp, _comm);
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		    (" elementToPrint  = ", &elementToPrint, statesStride) << "\n");

      if (_myRank == _ioRank) {
	for (CFuint i = 0; i < sendSize; ++i) {
	  if ((i+1)%statesStride > 0) {
	    fout->precision(16);
            fout->setf(ios::scientific,ios::floatfield);
	    *fout << elementToPrint[i] << " ";
	  }
	  else {
	    fout->precision(16);
            fout->setf(ios::scientific,ios::floatfield);
	    *fout << elementToPrint[i] << endl;
	  }
	}
      }

      //reset the all sendElement list to 0
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	sendElements[i] = 0;
	elementToPrint[i] = 0;
      }

      // update the count element for the current element type
      countElem += elementList.getSendDataSize(rangeID)/statesStride;
    }
  }

  CFLogInfo("States written \n");

  CFLogDebugMin( "ParCFmeshFileWriter::writeStateList() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileWriter::writeGeoList(CFuint iTRS, ofstream *const fout)
{
  if (_myRank == _ioRank) {
    *fout << "!LIST_GEOM_ENT" << "\n";
  }

  const vector<vector<CFuint> >&  trsInfo =
    MeshDataStack::getActive()->getTotalTRSInfo();

  SafePtr<vector<vector<vector<CFuint> > > > globalGeoIDS =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();

  const CFuint nbTRsInTRS = (*globalGeoIDS)[iTRS].size();
  CFLogDebugMin("nbTRsInTRS = " << nbTRsInTRS << "\n");

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
    CFLogDebugMin("nbGeosInLocalTR " << iTR << " is " << nbGeosInLocalTR << "\n");

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
		nbNodesStatesInTRGeo.size(), MPIStructDef::getMPIType(&nbNodesStatesInTRGeoTmp[0]), MPI_MAX, _comm);

  CFLogDebugMin("nbNodesStatesInTRGeo = " << nbNodesStatesInTRGeo << "\n");

  const CFuint nSend = _nbProc;
  const CFuint nbElementTypes = nbTRsInTRS;
  const CFuint nbLocalElements = trs->getLocalNbGeoEnts();

  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);

  // store global info about the global ID ranges for sending
  // and the size of each send

  CFuint maxElemSendSize = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint nbElementsInType = trsInfo[iTRS][iType];
    const CFuint maxNbNodesInType =  nbNodesStatesInTRGeo(iType, 0);
    const CFuint maxNbStatesInType = nbNodesStatesInTRGeo(iType, 1);
    // the max number of nodes and states includes two additional integers:
    // one for the number of nodes, the second for the number of states in
    // the current TR geometric entity
    const CFuint maxNodesPlusStatesData = maxNbNodesInType + maxNbStatesInType + 2;
    const CFuint maxElemSize = (nbElementsInType/nSend + nbElementsInType%nSend)*maxNodesPlusStatesData;
    const CFuint minElemSize = (nbElementsInType/nSend)*maxNodesPlusStatesData;

    // update the maximum possible element-list size to send
    maxElemSendSize = max(maxElemSendSize, maxElemSize);

    CFuint minElemID = 0;
    for (CFuint is = 0; is < nSend; ++is) {
      // first send per iType has maximum size
      // all the other sends (is != 0) per iType have minimum size
      const CFuint currSendSize = (is == 0) ? maxElemSize : minElemSize;

      // add the global size of the element data (each one of size=maxNodesPlusStates)
      // to communicate during this send
      elementList.addSendDataSize(currSendSize);
      const CFuint maxElemID = minElemID + currSendSize/maxNodesPlusStatesData;

      // insert the range
      elementList.addRange(minElemID, maxElemID);

      // set the new minimum as the old maximum
      minElemID = maxElemID;
    }
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
  
  CFuint rangeID = 0;
  for (CFuint iType = 0; iType < nbElementTypes; ++iType) {
    const CFuint maxNbNodesInType  = nbNodesStatesInTRGeo(iType, 0);
    const CFuint maxNbStatesInType = nbNodesStatesInTRGeo(iType, 1);
    const CFuint maxNodesPlusStates = maxNbNodesInType + maxNbStatesInType + 2;

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
	      CFLogDebugMin(_myRank << " maxNbNodesInType = " << maxNbNodesInType
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
	      CFLogDebugMin(_myRank << " state isend = " << isend << " , size = " << sendElements.size() << "\n");
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

      //   MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
      // 		 MPIStructDef::getMPIType(&sendElements[0]), MPI_MAX, _ioRank, _comm);
      
      MPI_Allreduce(&sendElements[0], &elementToPrint[0], maxElemSendSize, 
		    MPIStructDef::getMPIType(&elementToPrint[0]), MPI_MAX, _comm);
      
      CFLogDebugMax(_myRank << CFPrintContainer<vector<CFint> >
		    (" elementToPrint  = ", &elementToPrint, maxNodesPlusStates) << "\n");
      
      if (_myRank == _ioRank) {
	for (CFuint i = 0; i < sendSize; ++i) {
	  if ((i+1)%maxNodesPlusStates > 0) {
	    if (elementToPrint[i] != -1) {
	      *fout << elementToPrint[i] << " ";
	    }
	  }
	  else {
	    if (elementToPrint[i] != -1) {
	      *fout << elementToPrint[i] << endl;
	    }
	    else {
	      *fout << endl;
	    }
	  }
	}
      }

      //reset the all sendElement list to 0
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	sendElements[i] = -1;
	elementToPrint[i] = -1;
      }

      // update the count element for the current element type
      countElem += elementList.getSendDataSize(rangeID)/maxNodesPlusStates;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
