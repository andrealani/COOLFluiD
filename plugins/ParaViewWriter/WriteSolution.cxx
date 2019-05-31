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
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DataHandleOutput.hh"

#include "ParaViewWriter/ParaViewWriter.hh"
#include "ParaViewWriter/WriteSolution.hh"

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

MethodCommandProvider<WriteSolution, ParaWriterData, ParaViewWriterModule>
writeSolutionProvider("WriteSolution");

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write ParaView file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolution::WriteSolution(const std::string& name) : ParaWriterCom(name),
  socket_nodes("nodes"),
  socket_nstatesProxy("nstatesProxy")
{
  addConfigOptionsTo(this);

  m_fileFormatStr = "ASCII";
  setParameter("FileFormat",&m_fileFormatStr);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::execute()
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

const std::string WriteSolution::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::writeToBinaryFile()
{
 CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"Writing to binary files is not yet implemented in ParaViewWriter...");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "WriteSolution::writeToFileStream() => START\n");
  
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
  const CFuint nbrNodes = nodes.size();
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

  CFLog(VERBOSE, "WriteSolution::writeToFileStream() => nbVecComponents = "
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
  fout << "    <Piece NumberOfPoints=\"" << nbrNodes << "\" NumberOfCells=\"" << nbrCells << "\">\n";

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
      // get state in this node
      const RealVector& nodalState = *nodalStates.getState(iNode);

      // copy to temporary state variable
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        tempState[iEq] = nodalState[iEq];
      }

      // set temporary state ID
      const CFuint stateID = nodalStates.getStateLocalID(iNode);
      tempState.setLocalID(stateID);

      // set the node in the temporary state
      tempState.setSpaceCoordinates(nodes[iNode]);

      // dimensionalize the state
      updateVarSet->setDimensionalValues(tempState, dimState);

      // write the current variable to the file
      for (CFuint iVecComp = 0; iVecComp < nbVecComponents; ++iVecComp)
      {
        fout << scientific << setprecision(12) << dimState[vectorComponentIdxs[iVecComp]] << " ";
      }
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
        // get state in this node
        const RealVector& nodalState = *nodalStates.getState(iNode);

        // copy to temporary state variable
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          tempState[iEq] = nodalState[iEq];
        }

        // set temporary state ID
        const CFuint stateID = nodalStates.getStateLocalID(iNode);
        tempState.setLocalID(stateID);

        // set the node in the temporary state
        tempState.setSpaceCoordinates(nodes[iNode]);

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
      // get state in this node
      const RealVector& nodalState = *nodalStates.getState(iNode);

      // copy to temporary state variable
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        tempState[iEq] = nodalState[iEq];
      }

      // set temporary state ID
      const CFuint stateID = nodalStates.getStateLocalID(iNode);
      tempState.setLocalID(stateID);

      // set the node in the temporary state
      tempState.setSpaceCoordinates(nodes[iNode]);

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
        // get state in this node
        const RealVector& nodalState = *nodalStates.getState(iNode);

        // copy to temporary state variable
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          tempState[iEq] = nodalState[iEq];
        }

        // set temporary state ID
        const CFuint stateID = nodalStates.getStateLocalID(iNode);
        tempState.setLocalID(stateID);

        // set the node in the temporary state
        tempState.setSpaceCoordinates(nodes[iNode]);

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

      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
        fout << scientific << setprecision(12) << var(nodalStates.getStateLocalID(iNode), var_var, var_nbvars) << " ";

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
	    fout << scientific << setprecision(12) <<
	      var(iState, var_var, var_nbvars) << " ";
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
    for (CFuint iCoor = 0; iCoor < dim; ++iCoor)
    {
      fout << scientific << setprecision(12) << (*nodes[iNode])[iCoor]*refL << " ";
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
    const CFuint nbrNodes = cellNodes->nbCols(iCell);
    // node ordering for one cell is the same for VTK as in COOLFluiD
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      fout << (*cellNodes)(iCell,iNode) << " ";
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
    cellEndOffSet += cellNodes->nbCols(iCell);
    fout << cellEndOffSet << " ";
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
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // get VTK cell type
    const CFuint vtkCellType = getMethodData().getVTKCellTypeID((*elemType)[iElemType].getGeoShape(),(*elemType)[iElemType].getGeoOrder());

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
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

  CFLog(VERBOSE, "WriteSolution::writeToFileStream() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::writeBoundarySurface()
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

void WriteSolution::setup()
{
  CFAUTOTRACE;

  ParaWriterCom::setup();

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  updateVarSet->setup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolution::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_nstatesProxy);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD
