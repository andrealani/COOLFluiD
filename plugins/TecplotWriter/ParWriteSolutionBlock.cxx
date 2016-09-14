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
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/WriteListMap.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/ParWriteSolutionBlock.hh"

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

 MethodCommandProvider<ParWriteSolutionBlock, TecWriterData, TecplotWriterModule>
 parWriteSolutionBlockProvider("ParWriteSolutionBlock");

 //////////////////////////////////////////////////////////////////////////////
    
 ParWriteSolutionBlock::ParWriteSolutionBlock(const std::string& name) : 
   ParWriteSolution(name)
 {
   // change the default for m_onlyNodal
   // the user can still override this behaviour in the input file
   m_onlyNodal = false;
 }
    
 //////////////////////////////////////////////////////////////////////////////

 ParWriteSolutionBlock::~ParWriteSolutionBlock() 
 {
   CFAUTOTRACE;
 }

 //////////////////////////////////////////////////////////////////////////////

void ParWriteSolutionBlock::writeToBinaryFile()
{
 CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ParWriteSolutionBlock::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::setup() => start\n");
  
  ParWriteSolution::setup();
      
  CFLog(VERBOSE, "ParWriteSolutionBlock::setup() => end\n");
}      
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolutionBlock::writeNodeList(ofstream* fout, const CFuint iType,
					  SafePtr<TopologicalRegionSet> elements,
					  const bool isBoundary)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => start\n");
  
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  const bool printExtra = getMethodData().shouldPrintExtraValues();
  vector<string> extraVarNames = outputVarSet->getExtraVarNames();
  const CFuint nbExtraVars = extraVarNames.size();
  
  vector<bool> isVarNodal; 
  CFuint totNbVarsND = dim;
  for (CFuint i = 0; i < dim; ++i) {
    isVarNodal.push_back(true);
  }
  
  CFuint totNbVarsCC = 0;
  if (!getMethodData().onlyCoordinates()) {
    bool nodal = false;
    
    if (m_onlyNodal || isBoundary) {
      totNbVarsND += nbEqs;
      nodal = true;
    } 
    else if (!m_onlyNodal && !isBoundary) {
      totNbVarsCC += nbEqs;
      nodal = false;
    }
    
    for (CFuint i = 0; i < nbEqs; ++i) {
      isVarNodal.push_back(nodal);
    }
    
    if (printExtra) {
      if (m_onlyNodal || isBoundary) {
	totNbVarsND += nbExtraVars;
	nodal = true;
      }
      else if (!m_onlyNodal && !isBoundary) {
	totNbVarsCC += nbExtraVars;
	nodal = false;
      }
      
      for (CFuint i = 0; i < nbExtraVars; ++i) {
	isVarNodal.push_back(nodal);
      }
    }
    
    CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => 1 totNbVarsCC = " << totNbVarsCC << ", totNbVarsND = " << totNbVarsND << "\n");
    
    // data handles are ignored at the boundary
    if (!isBoundary) {
      // nodal data handle variables 
      totNbVarsND += datahandle_output->getVarNames().size();
      cf_assert(datahandle_output->getVarNames().size() == m_nodalvars.size());
      for (CFuint i = 0; i < datahandle_output->getVarNames().size(); ++i) {
	isVarNodal.push_back(true);
      }
      CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => 2 totNbVarsCC = " << totNbVarsCC << ", totNbVarsND = " << totNbVarsND << "\n");
      
      // cell-centered data handle variables 
      if (!m_onlyNodal) {
	totNbVarsCC += datahandle_output->getCCVarNames().size();
      }
      CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => 3 totNbVarsCC = " << totNbVarsCC << ", totNbVarsND = " << totNbVarsND << "\n");
      
      cf_assert(datahandle_output->getCCVarNames().size() == m_ccvars.size()); 
      for (CFuint i = 0; i < datahandle_output->getCCVarNames().size(); ++i) {
	isVarNodal.push_back(false);
      }
      CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => 4 totNbVarsCC = " << totNbVarsCC << ", totNbVarsND = " << totNbVarsND << "\n");
    }
  }
  cf_assert(totNbVarsND >= dim); // at least coordinates are printed in nodal mode
  cf_assert(isVarNodal.size() >= totNbVarsCC);
  cf_assert(isVarNodal.size() >= totNbVarsND);
  
  if (m_onlyNodal) {
    cf_assert(totNbVarsCC == 0);
  }
  else {
    cf_assert(isVarNodal.size() == totNbVarsCC + totNbVarsND);
  }
  
  const CFuint nodesStride = 1; // one variable at a time
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  cf_assert(states.size() > 0);
  
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();
  
  // this is a sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];
  
  TeclotTRSType& tt = *_mapTrsName2TecplotData.find(elements->getName());
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = 1; // just nodes
  
  RealVector dimState(nbEqs);
  // size should be set in the VarSet but for safety we resize it here too
  RealVector extraValues;
  if (nbExtraVars > 0) extraValues.resize(nbExtraVars);
  State tempState;
  
  const string writerName = getMethodData().getNamespace() + "_Writers";
  Group& wg = PE::GetPE().getGroup(writerName);
  
  // fill in the writer list for nodal variables
  const CFuint totNbNodes = tt.totalNbNodesInType[iType];
  const CFuint nbLocalElementsND = tt.nodesInType[iType].size();
  CFuint totalToSendND = 0;
  WriteListMap elementListND;
  elementListND.reserve(nbElementTypes, nSend, nbLocalElementsND);
  elementListND.fill(totNbNodes, nodesStride, totalToSendND);
  
  SafePtr<vector<ElementTypeData> > me =
    MeshDataStack::getActive()->getElementTypeData();
  
  // fill in the writer list for cell-centered variables
  const CFuint nbElemsInType = (*me)[iType].getNbTotalElems();
  CFuint nbLocalElementsCC = 0;
  CFuint totalToSendCC = 0;
  WriteListMap elementListCC;
  if (!m_onlyNodal && !isBoundary) { 
    SafePtr<vector<ElementTypeData> > elementType =
      MeshDataStack::getActive()->getElementTypeData(elements->getName());
    ElementTypeData& eType = (*elementType)[iType];
    
    nbLocalElementsCC = eType.getNbElems();
    elementListCC.reserve(nbElementTypes, nSend, nbLocalElementsCC);
    elementListCC.fill(nbElemsInType, nodesStride, totalToSendCC);
  }
  
  const CFuint wordFormatSize = 22;
  
  vector<MPI_Offset> wOffset(_nbWriters, tt.headerOffset[iType][1]); 
  // set the offsets for the nodes of this element type
  tt.nodesOffset[iType].first  = tt.headerOffset[iType][1];
  const CFuint totalToSend = totNbVarsND*totalToSendND + totNbVarsCC*totalToSendCC;
  tt.nodesOffset[iType].second = tt.nodesOffset[iType].first + totalToSend*wordFormatSize;
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => offsets = [" 
	<<  tt.nodesOffset[iType].first << ", " << tt.nodesOffset[iType].second << "]\n");
  
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = std::max(elementListND.getMaxElemSize(),
					  elementListCC.getMaxElemSize());
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  const vector<CFuint>& nodesInType = tt.nodesInType[iType];
  for (CFuint iElem = 0; iElem < nbLocalElementsND; ++iElem) {
    const CFuint globalTypeID = tt.mapNodeID2NodeIDByEType[iType]->find(nodesInType[iElem]);
    elementListND.insertElemLocalID(iElem, globalTypeID, 0);
  }
  elementListND.endElemInsertion(_myRank);
  
  SafePtr< vector<CFuint> > globalElementIDs = MeshDataStack::getActive()->getGlobalElementIDs();
  
  if (!m_onlyNodal && !isBoundary) {
    SafePtr<vector<ElementTypeData> > elementType =
      MeshDataStack::getActive()->getElementTypeData(elements->getName());
    ElementTypeData& eType = (*elementType)[iType];
    // insert in the write list the local IDs of the elements
    // the range ID is automatically determined inside the WriteListMap
    CFuint elemID = eType.getStartIdx();
    for (CFuint iElem = 0; iElem < nbLocalElementsCC; ++iElem, ++elemID) {
      elementListCC.insertElemLocalID(elemID, (*globalElementIDs)[elemID], 0);
    }
    elementListCC.endElemInsertion(_myRank);
  }
  
  // buffer data to send
  vector<CFreal> sendElements(maxElemSendSize, 0);
  vector<CFreal> elementToPrint(maxElemSendSize, 0);
  
  const CFuint totNbVars = totNbVarsND + totNbVarsCC;
  if (getMethodData().onlyCoordinates()) {
    cf_always_assert(totNbVarsND == dim);
    cf_always_assert(totNbVars == dim); 
  }
  
  const CFuint nbDHNDVars  = datahandle_output->getVarNames().size();
  const CFuint nbDHCCVars  = datahandle_output->getCCVarNames().size();
  const CFuint endNbEqs = dim + nbEqs;
  const CFuint endNbExtraVars = (printExtra) ? endNbEqs + nbExtraVars : endNbEqs;
  const CFuint endNbDHNDVars = endNbExtraVars + nbDHNDVars;
  const CFuint endNbDHCCVars = endNbDHNDVars + nbDHCCVars;
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => endNbEqs       = " << endNbEqs << "\n"); 
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => endNbExtraVars = " << endNbExtraVars << "\n"); 
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => endNbDHNDVars  = " << endNbDHNDVars << "\n");
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => endNbDHCCVars  = " << endNbDHCCVars << "\n");
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => totNbVars      = " << totNbVars << "\n"); 
  
  if (!m_onlyNodal) {
    cf_assert((!isBoundary && endNbDHCCVars == totNbVars) || 
	      (isBoundary && endNbExtraVars == totNbVars));
  }
  
  for (CFuint iVar = 0; iVar < totNbVars; ++iVar) {
    CFint wRank = -1; 
    CFuint wSendSize = 0;
    CFuint rangeID = 0;
    CFuint countElem = 0;
    
    cf_assert(iVar < isVarNodal.size());
    const bool isNodal = isVarNodal[iVar]; 
    CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => isNodal[" << iVar << "] = " << isNodal << " \n"); 
    
    // if the output has to be NODAL, discard CC variables
    if (!m_onlyNodal || (m_onlyNodal && isNodal)) {
      WriteListMap& elementList = (isNodal) ? elementListND : elementListCC;
      
      for (CFuint is = 0; is < nSend; ++is, ++rangeID) {
	bool isRangeFound = false;
	WriteListMap::List elist = elementList.find(rangeID, isRangeFound);
	
	if (isRangeFound) {
	  CFuint eSize = 0;
	  for (WriteListMap::ListIterator it = elist.first; it != elist.second; ++it, ++eSize) {
	    const CFuint localElemID = it->second;
	    CFuint globalElemID = 0;
	    CFuint dofID = 0;
	    bool isUpdatable = false;
	    if (isNodal) {
	      globalElemID = tt.mapNodeID2NodeIDByEType[iType]->find(nodesInType[localElemID]);
	      dofID = _mapGlobal2LocalNodeID.find(nodesInType[localElemID]);
	      if (nodes[dofID]->isParUpdatable()) isUpdatable = true;
	    }
	    else {
	      globalElemID = (*globalElementIDs)[localElemID];
	      dofID = localElemID; // AL: double check this!!!
	      if (states[dofID]->isParUpdatable()) isUpdatable = true;
	    }
	    
	    const CFuint sendElemID = globalElemID - countElem;
	    
	    // this fix has to be added EVERYWHERE when writing states in parallel
	    if (isUpdatable) {
	      CFuint isend = sendElemID*nodesStride;
	      
	      if (iVar < dim) {
		cf_assert(isend < sendElements.size());
		cf_assert(dofID < nodes.size());
		sendElements[isend++] = (*nodes[dofID])[iVar]*refL;
		// break;
	      }
	      
	      if (isNodal && iVar >= dim) { 
		const CFuint stateID = nodalStates.getStateLocalID(dofID);
		tempState.setLocalID(stateID);
		// the node is set  in the temporary state
		tempState.setSpaceCoordinates(nodes[dofID]);
		
		if (iVar < endNbExtraVars) {
		  const RealVector& currState = *nodalStates.getState(dofID);
		  for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
		    tempState[ieq] = currState[ieq];
		  }
		  
		  if (iVar >= dim && iVar < endNbEqs) {
		    cf_assert(isend < sendElements.size());
		    outputVarSet->setDimensionalValues(tempState, dimState);
		    sendElements[isend++] = dimState[iVar-dim];
		  }
		  
		  if (printExtra && (iVar >= endNbEqs && iVar < endNbExtraVars)) {
		    // this is EXTREMELY inefficient !!!! 
		    // setDimensionalValuesPlusExtraValues() will be called nbExtraVars times per state !!!
		    // dimensionalize the solution
		    outputVarSet->setDimensionalValuesPlusExtraValues(tempState, dimState, extraValues);
		    cf_assert(isend < sendElements.size());
		    sendElements[isend++] = extraValues[iVar-endNbEqs];
		  }
		}
		
		if (!isBoundary) {
		  if (nbDHNDVars > 0 && iVar >= endNbExtraVars && iVar < endNbDHNDVars) {
		    datahandle_output->fillStateData(&sendElements[0], stateID, isend, iVar-endNbExtraVars);
		  }
		}
	      }
	      else if ((!isNodal) && iVar >= dim) {
		// cell centered case
		if (iVar >= dim && iVar < endNbEqs) {
		  const State& currState = *states[dofID];
		  outputVarSet->setDimensionalValues(currState, dimState);
		  cf_assert(isend < sendElements.size());
		  sendElements[isend++] = dimState[iVar-dim];
		}
		
		if (printExtra && (iVar >= endNbEqs && iVar < endNbExtraVars)) {
		  // this is EXTREMELY inefficient !!!! 
		  // setDimensionalValuesPlusExtraValues() will be called nbExtraVars times per state !!!
		  // dimensionalize the solution
		  const State& currState = *states[dofID];
		  outputVarSet->setDimensionalValuesPlusExtraValues(currState, dimState, extraValues);
		  cf_assert(isend < sendElements.size());
		  cf_assert(iVar-endNbEqs < extraValues.size());
		  sendElements[isend++] = extraValues[iVar-endNbEqs]; 
		}
		
		if (nbDHCCVars > 0 && iVar >= endNbDHNDVars) {
		  datahandle_output->fillStateDataCC(&sendElements[0], dofID, isend, iVar - endNbDHNDVars);
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
	CFLog(DEBUG_MIN, "ParWriteSolutionBlock::writeNodeList() => P[" << _myRank << "] => offset = " << wOffset[wRank] << "\n");
	// each writer can now concurrently write all the collected data (related to one element type)
	cf_assert(wRank >= 0);
	
	// point to the corresponding writing location
	fout->seekp(wOffset[wRank]);
	CFLog(VERBOSE, "For iVar=" << iVar << " seeking position wOffset[wRank] = " << wOffset[wRank] << "\n");
	CFLog(VERBOSE, "wSendSize = " << wSendSize << ", nodesStride = " << nodesStride << "\n");
	// const CFuint sizeLine = 8;
	for (CFuint i = 0; i < wSendSize; ++i) {
	  // this fix is needed for ensuring consistency in the format and avoid 
	  // to have entries with 3 digits in the exponent (e.g. +A.BC..e-XXX)
	  if (std::abs(elementToPrint[i]) < 1e-50) elementToPrint[i] = 0.;
	  
	  //	const CFreal coeff = (countN < dim) ? refL : 1.;
	  // this format corresponds to a line of 22*nodesStride bytes  
	  fout->precision(14);
	  fout->setf(ios::scientific,ios::floatfield); 
	  fout->setf(ios::showpos);
	  // if ((i+1)%sizeLine > 0) {
	  // *fout << elementToPrint[i] << " ";
	  //	}
	  //else {
	  *fout << elementToPrint[i] << "\n"; 
	  // }
	}
	
	MPI_Offset lastpos = fout->tellp();
	MPI_Offset maxpos  = 0;
	MPI_Allreduce(&lastpos, &maxpos, 1, MPIStructDef::getMPIOffsetType(), MPI_MAX, wg.comm);
	
	if (iVar == totNbVars-1) {
	  if (tt.nodesOffset[iType].second != maxpos) {
	    CFLog(ERROR, "tt.nodesOffset[iType].second (" << tt.nodesOffset[iType].second 
		  << ") != maxpos (" << maxpos << ")" << "\n");
	    cf_assert(tt.nodesOffset[iType].second == maxpos);
	  }
	}
	else {
	  // update the starting offset
	  wOffset.assign(_nbWriters, maxpos); 
	}
	fout->seekp(maxpos);
      }
      
      //reset the all sendElement list to 0
      for (CFuint i = 0; i < maxElemSendSize; ++i) {
	elementToPrint[i] = 0;
      }
    }
  }
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeNodeList() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolutionBlock::writeZoneHeader(std::ofstream* fout, 
					    const CFuint iType,
					    const string& geoShape,
					    const CFuint nbNodesInType,
					    const CFuint nbElemsInType,
					    const string& geoType,
					    const string& end,
					    const bool isBoundary) 
{
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeZoneHeader() => start\n");
  
  *fout << "ZONE "
	<< "  T= \"ZONE" << iType << " " << geoShape <<"\""
	<< ", N=" << std::noshowpos << nbNodesInType
	<< ", E=" << std::noshowpos << nbElemsInType
        << ", ZONETYPE=FE" << geoType
	<< ", DATAPACKING=BLOCK, STRANDID=1, SOLUTIONTIME=";
  fout->precision(14);  
  fout->setf(ios::scientific,ios::floatfield);
  const CFreal timeDim = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  const CFreal nbIter  = (CFreal)SubSystemStatusStack::getActive()->getNbIter();
  *fout << ((timeDim > 0.) ? timeDim : nbIter);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  *fout << ", VARLOCATION=( ";
  if (dim!=1) {
    *fout   << "[" << 1 << "-" << dim << "]=NODAL";
  }
  else {
    *fout << dim << "=NODAL";
  }
  
  if (!getMethodData().onlyCoordinates()) {
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    string locationString = "NODAL";
    if (!isBoundary) {
      locationString = m_onlyNodal ? "NODAL" : "CELLCENTERED";
    }
    
    if (nbEqs!=1) {
      *fout << ",[" << dim+1 << "-" << dim+nbEqs << "]=" << locationString;
    }
    else {
      *fout << "," << dim+1 << "=" << locationString;
    }
    
    CFuint nbExtraVars = 0;
    if (getMethodData().shouldPrintExtraValues()) { 
      SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
      vector<string> extraVarNames = outputVarSet->getExtraVarNames();
      nbExtraVars = extraVarNames.size();
      if (nbExtraVars > 0) {
	*fout << ",[" << dim+nbEqs+1 << "-" << dim+nbEqs+nbExtraVars << "]=" << locationString;
      }
      else {
	*fout << "," << dim+nbEqs+1 << "=" << locationString;
      }
    }
    
    if (!isBoundary) {
      if (!m_nodalvars.empty()) {
	if (m_nodalvars.size()!=1) {
	  *fout << ",[" << dim+nbEqs+1+nbExtraVars << "-" << 
	    dim+nbEqs+nbExtraVars+m_nodalvars.size() << "]=NODAL";
	}
	else {
	  *fout << "," << dim+nbEqs+1+nbExtraVars << "=NODAL";
	}
      }
      
      // only consider DataHandle CC if !m_onlyNodal
      if (!m_ccvars.empty() && !m_onlyNodal) {
	if(m_ccvars.size()!=1) {
	  *fout << ",[" << dim+nbEqs+nbExtraVars+m_nodalvars.size()+1 << "-" << 
	    dim+nbEqs+nbExtraVars+m_nodalvars.size()+m_ccvars.size() << "]=CELLCENTERED";
	}
	else {
	  *fout << "," << dim+nbEqs+nbExtraVars+m_nodalvars.size()+1 << "=CELLCENTERED";
	}
      }
    }
  }
  
  *fout << ")" << end;
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeZoneHeader() => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////
  
  } // namespace TecplotWriter
  
} // namespace COOLFluiD
