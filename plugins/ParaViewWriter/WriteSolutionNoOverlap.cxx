// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include <iomanip>

#include "Common/CFMap.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"

#include "ParaViewWriter/ParaViewWriter.hh"
#include "ParaViewWriter/WriteSolutionNoOverlap.hh"

#include "Common/OSystem.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteSolutionNoOverlap, ParaWriterData, ParaViewWriterModule>
writeSolutionNoOverlapProvider("WriteSolutionNoOverlap");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionNoOverlap::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write ParaView file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionNoOverlap::WriteSolutionNoOverlap(const std::string& name) :
  ParaWriterCom(name),
  socket_nodes("nodes"),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states"),
  socket_gstates("gstates"),
  m_nbCellsNoPartition(0),
  m_isPartitionNode(),
  m_cellWithPartitionNodes(),
  m_oldToNewNodeID(),
  m_newToOldNodeID(),
  m_nbCellsNoPartitionInElemType()
{
  addConfigOptionsTo(this);

  m_fileFormatStr = "ASCII";
  setParameter("FileFormat",&m_fileFormatStr);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionNoOverlap::execute()
{
  CFLog(INFO, "Writing solution to: " << getMethodData().getFilename().string() << "\n");
  
  if(m_fileFormatStr == "ASCII")
  {
    writeToFile(getMethodData().getFilename());
  }
  else
  {
    cf_assert(m_fileFormatStr == "BINARY");

    writeToBinaryFile();
  }

}

//////////////////////////////////////////////////////////////////////////////

const std::string WriteSolutionNoOverlap::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionNoOverlap::writeToBinaryFile()
{
 CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"Writing to binary files is not yet implemented in ParaViewWriter...");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionNoOverlap::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "WriteSolutionNoOverlap::writeToFileStream() => START\n");
  
  if (!getMethodData().onlySurface())
  {
    // get the nodes datahandle
    DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
    
    // get iterator for nodal states datahandle
    DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy = socket_nstatesProxy.getDataHandle();

    // this is a sort of handle for the nodal states
    // (which can be stored as arrays of State*, RealVector* or
    // RealVector but they are used as arrays of RealVector*)
    ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];
      
  // get the cells
  SafePtr<TopologicalRegionSet> elements = MeshDataStack::getActive()->getTrs("InnerCells");

  // get cell-node connectivity
  SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // number of cells and nodes
  const CFuint nbrCells = elements->getLocalNbGeoEnts();
  const CFuint nbrNodes = m_newToOldNodeID.size();
  cf_assert(nbrNodes <= nodes.size());
  cf_assert(cellNodes->nbRows() == nbrCells);
  
  // we will assume that the number of nodes is the same as
  // the number of states but the connectivity might be different
  //   cf_assert(nodes.size() == nodalStates.getSize());
  
  // get the element type data
  SafePtr<vector<ElementTypeData> > elemType =  MeshDataStack::getActive()->getElementTypeData();

  // get dimensionality, number of variables and reference length
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  
  // variable that holds the indices of the variables that are vector components
  vector<CFuint> vectorComponentIdxs;
  /// @note this is a rather ugly piece of code, a check is made on the number of equations
  /// to avoid that vector components are searched when the physical model is a linear advection for instance
  /// this piece of code puts the velocity components (or the momentum components) in a vector
  /// the magnetic inductance vector B in the case of MHD is not put in a vector here like this!
  /// AL: fix to get the velocity IDs directly from the physics
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  updateVarSet->setStateVelocityIDs(vectorComponentIdxs);
  
  // variable tht holds the indices of the scalar variables
  const CFuint nbVecComponents = vectorComponentIdxs.size();
  const CFuint nbScalars = nbEqs-nbVecComponents;

  CFLog(VERBOSE, "WriteSolutionNoOverlap::writeToFileStream() => nbVecComponents = "
	<< nbVecComponents << ",  nbScalars = " << nbScalars << "\n");
  
  vector<CFuint> scalarVarIdxs(nbScalars);
  CFuint iScalar = 0;
  for  (CFuint iEq = 0; iEq < nbEqs; ++iEq)
  {
    bool addIdx = true;
    for (CFuint iVecComp = 0; iVecComp < nbVecComponents; ++iVecComp)
    {
      if (iEq == vectorComponentIdxs[iVecComp])
      {
        addIdx = false;
      }
    }
    if (addIdx)
    {
      scalarVarIdxs[iScalar] = iEq;
      ++iScalar;
    }
  }
  cf_assert(iScalar == nbScalars);
  
  // get variable names
  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  // open VTKFile element
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=";
  if (isLittleEndian())
  {
    fout << "\"LittleEndian\">\n";
  }
  else
  {
    fout << "\"BigEndian\">\n";
  }

  // open UnstructuredGrid element
  fout << "  <UnstructuredGrid>\n";

  // open Piece element
  fout << "    <Piece NumberOfPoints=\"" << nbrNodes << "\" NumberOfCells=\"" << m_nbCellsNoPartition << "\">\n";

  // open PointData element
//   fout << "      <PointData>\n";
  fout << "   <PointData Scalars=\"" << varNames[0] << "\">\n";

  // some helper states
  RealVector dimState(nbEqs);
  RealVector extraValues; // size will be set in the VarSet
  State tempState;

  // write the (velocity or momentum) vectors
  if ((nbVecComponents > 0) && (!getMethodData().writeVectorAsComponents()))
  {
    cf_assert(nbVecComponents >= 2);
    // open DataArray element
    fout << "        <DataArray Name=\"" << varNames[vectorComponentIdxs[1]] << "\" NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";
    fout << "          ";

    // loop over nodes
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      const CFuint nodeID = m_newToOldNodeID[iNode];
      
      // get state in this node
      const RealVector& nodalState = *nodalStates.getState(nodeID);

      // copy to temporary state variable
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        tempState[iEq] = nodalState[iEq];
      }
      
      // set temporary state ID
      const CFuint stateID = nodalStates.getStateLocalID(nodeID);
      tempState.setLocalID(stateID);
      
      // set the node in the temporary state
      tempState.setSpaceCoordinates(nodes[nodeID]);
      
      // dimensionalize the state
      updateVarSet->setDimensionalValues(tempState, dimState);

      // write the current variable to the file
      for (CFuint iVecComp = 0; iVecComp < nbVecComponents; ++iVecComp)
      {
        fout << scientific << setprecision(12) << dimState[vectorComponentIdxs[iVecComp]] << " ";
	// cout << "[" << vectorComponentIdxs[iVecComp] << "] " << scientific << setprecision(12) << dimState[vectorComponentIdxs[iVecComp]] << " ";
      }
      // cout << endl;
      for (CFuint iVecComp = nbVecComponents; iVecComp < 3; ++iVecComp)
      {
	fout << scientific << setprecision(1) << 0.0 << " ";
      }

    }
    fout << "\n";

    // close DataArray element
    fout << "        </DataArray>\n";
  } else if (nbVecComponents > 0) {

    for (CFuint iVecComp = 0; iVecComp < nbVecComponents; ++iVecComp)
    {

      fout << "        <DataArray Name=\"" << varNames[vectorComponentIdxs[iVecComp]] << "\" type=\"Float32\" format=\"ascii\">\n";
      fout << "          ";

      // loop over nodes
      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
      {
	const CFuint nodeID = m_newToOldNodeID[iNode];
	
        // get state in this node
        const RealVector& nodalState = *nodalStates.getState(nodeID);

        // copy to temporary state variable
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          tempState[iEq] = nodalState[iEq];
        }

        // set temporary state ID
        const CFuint stateID = nodalStates.getStateLocalID(nodeID);
        tempState.setLocalID(stateID);

        // set the node in the temporary state
        tempState.setSpaceCoordinates(nodes[nodeID]);

        // dimensionalize the state
        updateVarSet->setDimensionalValues(tempState, dimState);

        // write the current variable to the file
        fout << scientific << setprecision(12) << dimState[vectorComponentIdxs[iVecComp]] << " ";

      }
      fout << "\n";

      // close DataArray element
      fout << "        </DataArray>\n";
    }
  }

  // loop over the scalars
  for (CFuint iScalar = 0; iScalar < nbScalars; ++iScalar)
  {
    // index of this scalar'
    const CFuint iVar = scalarVarIdxs[iScalar];

    // open DataArray element
    fout << "        <DataArray Name=\"" << varNames[iVar] << "\" type=\"Float32\" format=\"ascii\">\n";
    fout << "          ";

    // loop over nodes
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      const CFuint nodeID = m_newToOldNodeID[iNode];
      
      // get state in this node
      const RealVector& nodalState = *nodalStates.getState(nodeID);

      // copy to temporary state variable
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        tempState[iEq] = nodalState[iEq];
      }

      // set temporary state ID
      const CFuint stateID = nodalStates.getStateLocalID(nodeID);
      tempState.setLocalID(stateID);

      // set the node in the temporary state
      tempState.setSpaceCoordinates(nodes[nodeID]);

      // dimensionalize the state
      updateVarSet->setDimensionalValues(tempState, dimState);

      // write the current variable to the file
      fout << scientific << setprecision(12) << dimState[iVar] << " ";
    }
    fout << "\n";

    // close DataArray element
    fout << "        </DataArray>\n";
  }

  // if extra variables are to be outputted
  if (getMethodData().printExtraValues())
  {
    // variable for extra values
    RealVector extraValues; // size will be set in the VarSet

    // extra variable names
    const vector<std::string>& extraVarNames = updateVarSet->getExtraVarNames();
    const CFuint nbrExtraVars = extraVarNames.size();

    // loop over the extra variables
    for (CFuint iVar = 0 ;  iVar < nbrExtraVars; ++iVar)
    {
      // open DataArray element
      fout << "        <DataArray Name=\"" << extraVarNames[iVar] << "\" type=\"Float32\" format=\"ascii\">\n";
      fout << "          ";

      // loop over nodes
      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
      {
	const CFuint nodeID = m_newToOldNodeID[iNode];
        // get state in this node
        const RealVector& nodalState = *nodalStates.getState(nodeID);

        // copy to temporary state variable
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          tempState[iEq] = nodalState[iEq];
        }

        // set temporary state ID
        const CFuint stateID = nodalStates.getStateLocalID(nodeID);
        tempState.setLocalID(stateID);

        // set the node in the temporary state
        tempState.setSpaceCoordinates(nodes[nodeID]);

        // dimensionalize the state
        updateVarSet->setDimensionalValuesPlusExtraValues(tempState, dimState, extraValues);

        // write the current variable to the file
        fout << scientific << setprecision(12) << extraValues[iVar] << " ";
      }
      fout << "\n";

      // close DataArray element
      fout << "        </DataArray>\n";
    }
  }


  // print datahandles with state based data
  {
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    datahandle_output->getDataHandles();
    std::vector< std::string > dh_varnames = datahandle_output->getVarNames();
    // loop over the state based variables

    for (CFuint iVar = 0; iVar < dh_varnames.size(); ++iVar)
    {
      // open DataArray element
      fout << "        <DataArray Name=\"" << dh_varnames[iVar] << "\" type=\"Float32\" format=\"ascii\">\n";
      fout << "          ";

      DataHandleOutput::DataHandleInfo var_info = datahandle_output->getStateData(iVar);
      CFuint var_var = var_info.first;
      CFuint var_nbvars = var_info.second;
      DataHandle<CFreal> var = var_info.third;

      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode) {
	const CFuint nodeID = m_newToOldNodeID[iNode];
	fout << scientific << setprecision(12) << var(nodalStates.getStateLocalID(nodeID), var_var, var_nbvars) << " ";
      }
      fout << "\n";

      // close DataArray element
      fout << "        </DataArray>\n";
    }
  }

  // close PointData element
  fout << "      </PointData>\n";

  // print datahandles with state based data
  {
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    datahandle_output->getDataHandles();
    std::vector< std::string > dh_varnames = datahandle_output->getCCVarNames();
    // loop over the state based variables
    
    if (dh_varnames.size() > 0) {
      // cell-based data
      // open CellData element
      fout << "   <CellData Scalars=\"" << dh_varnames[0] << "\">\n";
      
      for (CFuint iVar = 0; iVar < dh_varnames.size(); ++iVar)
	{
	  // open DataArray element
	  fout << "        <DataArray Name=\"" << dh_varnames[iVar] << "\" type=\"Float32\" format=\"ascii\">\n";
	  fout << "          ";
	  
	  DataHandleOutput::DataHandleInfo var_info = datahandle_output->getCCData(iVar);
	  CFuint var_var = var_info.first;
	  CFuint var_nbvars = var_info.second;
	  DataHandle<CFreal> var = var_info.third;
	  
	  for (CFuint iState = 0; iState < nbrCells; ++iState) {
          if (!m_cellWithPartitionNodes[iState]) {
	    fout << scientific << setprecision(12) << var(iState, var_var, var_nbvars) << " ";
	   }
          }	  
	  fout << "\n";
	  
	  // close DataArray element
	  fout << "        </DataArray>\n";
	}
      
      // close CellData element
      fout << "      </CellData>\n";
    }
  }
  
  // open Points element
  fout << "      <Points>\n";

  // open DataArray element
  fout << "        <DataArray NumberOfComponents=\"" << 3 << "\" type=\"Float32\" format=\"ascii\">\n";
  fout << "          ";

  // loop over nodes to write coordinates
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    const CFuint nodeID = m_newToOldNodeID[iNode];
    for (CFuint iCoor = 0; iCoor < dim; ++iCoor)
    {
      fout << scientific << setprecision(12) << (*nodes[nodeID])[iCoor]*refL << " ";
    }
    for (CFuint iCoor = dim; iCoor < 3; ++iCoor)
    {
      fout << scientific << setprecision(1) << 0.0 << " ";
    }
  }
  fout << "\n";

  // close DataArray element
  fout << "        </DataArray>\n";

  // close Points element
  fout << "      </Points>\n";

  // open Cells element
  fout << "      <Cells>\n";

  // open DataArray element (for cell-node connectivity)
  fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  fout << "          ";

  // loop over element types to write cell-node connectivity
  for (CFuint iCell = 0; iCell < nbrCells; ++iCell)
  {
    if (!m_cellWithPartitionNodes[iCell]) {
    const CFuint nbrNodes = cellNodes->nbCols(iCell);
    // node ordering for one cell is the same for VTK as in COOLFluiD
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      const CFuint newNodeID = m_oldToNewNodeID[(*cellNodes)(iCell,iNode)];
      fout << newNodeID << " ";
    }
   }
  }
  fout << "\n";

  // close DataArray element (for cell-node connectivity)
  fout << "        </DataArray>\n";

  // open DataArray element (for offsets in cell-node connectivity)
  fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  fout << "          ";

  // loop over element types to write offsets in cell-node connectivity (offset of the end of the connectivity for each cell)
  CFuint cellEndOffSet = 0;
  for (CFuint iCell = 0; iCell < nbrCells; ++iCell)
  {
    if (!m_cellWithPartitionNodes[iCell]) {
     cellEndOffSet += cellNodes->nbCols(iCell);
     fout << cellEndOffSet << " ";
    }
  }
  fout << "\n";

  // close DataArray element (for offsets in cell-node connectivity)
  fout << "        </DataArray>\n";

  // open DataArray element (for cell types)
  fout << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  fout << "          ";

  // loop over element types to write cell types
  /// @warning (element indexes (elemIdx) should increase monotonically here in order for this to be correct!!!)
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get VTK cell type
    const CFuint vtkCellType = getMethodData().getVTKCellTypeID
      ((*elemType)[iElemType].getGeoShape(),(*elemType)[iElemType].getGeoOrder());

    // loop over cells
    const CFuint nbCellsNoPartitionInElemType = m_nbCellsNoPartitionInElemType[iElemType];
    for (CFuint elemIdx = 0; elemIdx < nbCellsNoPartitionInElemType; ++elemIdx) {
      fout << vtkCellType << " ";
    }
  }
  fout << "\n";

  // close DataArray element (for cell types)
  fout << "        </DataArray>\n";

  // close Cells element
  fout << "      </Cells>\n";

  // close Piece element
  fout << "    </Piece >\n";

  // close UnstructuredGrid element
  fout << "  </UnstructuredGrid>\n";

  // close VTKFile element
  fout << "</VTKFile>\n";

  // close the file
  fout.close();

  } // if only surface

  // write boundary surface data
  writeBoundarySurface();

  CFLog(VERBOSE, "WriteSolutionNoOverlap::writeToFileStream() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionNoOverlap::writeBoundarySurface()
{
  CFAUTOTRACE;

  // check if there are surface TRSs to write
  if (getMethodData().getSurfaceTRSsToWrite().empty())
    return;

  // get surface TRS list
  const std::vector<std::string>& surfTRS = getMethodData().getSurfaceTRSsToWrite();

  // count number of surface TRs to write
  CFuint countTRToWrite = 0;
  std::vector<std::string>::const_iterator itr = surfTRS.begin();
  for(; itr != surfTRS.end(); ++itr)
  {
    Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);
    const CFuint nbTRs = currTrs->getNbTRs();
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
    {
      SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
      if (tr->getLocalNbGeoEnts() > 0)
      {
        ++countTRToWrite;
      }
    }
  }

  CFLog(WARN,"Writing of boundary surface data is not yet implemented in ParaViewWriter...");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionNoOverlap::setup()
{
  CFAUTOTRACE;

  ParaWriterCom::setup();

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  updateVarSet->setup();

  // get the nodes datahandle
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  const CFuint nbAllNodes = nodes.size();
  m_isPartitionNode.resize(nbAllNodes, false);
  
  GeometricEntityPool<FaceTrsGeoBuilder> faceBuilder;
  faceBuilder.setup();
  faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = faceBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> pTRS = MeshDataStack::getActive()->getTrs("PartitionFaces");
  faceData.trs = pTRS;
  faceData.isBFace = true;
  
  const CFuint nbFaces = faceData.trs->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    faceData.idx = iFace;
    GeometricEntity *const face = faceBuilder.buildGE();
    const CFuint nbNodes = face->nbNodes();
    const vector<Node*>& nodes = *face->getNodes();
    for (CFuint in = 0; in < nbNodes; ++in) {
      const CFuint nodeID = nodes[in]->getLocalID();
      m_isPartitionNode[nodeID] = true;
    }
    faceBuilder.releaseGE();
  }
  
  m_newToOldNodeID.reserve(nbAllNodes);
  m_oldToNewNodeID.resize(nbAllNodes, -1);
  
  CFuint newNodeID = 0;
  for (CFuint iNode = 0; iNode < nbAllNodes; ++iNode) {
    if (!m_isPartitionNode[iNode]) {
      cf_assert(newNodeID <= iNode);
      m_oldToNewNodeID[iNode] = newNodeID;
      m_newToOldNodeID.push_back(iNode);
      newNodeID++;
      cf_assert(m_newToOldNodeID.size() == newNodeID);
    }
    else {
      m_oldToNewNodeID[iNode] = -1;
    }
  }
  
  cf_assert(m_newToOldNodeID.size() <= nbAllNodes);
  
  // get cell-node connectivity
  SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  const CFuint nbAllCells = cellNodes->nbRows();
  m_cellWithPartitionNodes.resize(nbAllCells, false);

  m_nbCellsNoPartition = 0;

  // get the element type data
  SafePtr<vector<ElementTypeData> > elemType =  MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();
  m_nbCellsNoPartitionInElemType.resize(nbrElemTypes, 0);
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType) {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx) {
      const CFuint nbrNodes = cellNodes->nbCols(elemIdx);
      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode) {
	const CFuint nodeID = (*cellNodes)(elemIdx,iNode);
	if (m_isPartitionNode[nodeID]) {
	  m_cellWithPartitionNodes[elemIdx] = true;
	  break;
	}
      }
      if (!m_cellWithPartitionNodes[elemIdx]) {
	m_nbCellsNoPartition++;
	m_nbCellsNoPartitionInElemType[iElemType]++;
      }
    }
  }
  
  // AL: loop over element types to write cell-node connectivity for hybrid meshes
  /*for (CFuint iCell = 0; iCell < nbAllCells; ++iCell) {
    const CFuint nbrNodes = cellNodes->nbCols(iCell);
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode) {
1      const CFuint nodeID = (*cellNodes)(iCell,iNode);
      if (m_isPartitionNode[nodeID]) {
	m_cellWithPartitionNodes[iCell] = true;
	m_nbCellsNoPartition++;
	break;
      }
    }
    } */
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionNoOverlap::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD
