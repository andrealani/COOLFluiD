// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include <iomanip>

#include "Common/COOLFluiD.hh"

#include "pcgnslib.h"
//#include <sys/stat.h>

#include "Common/PE.hh"
#include "Common/MPI/MPIIOFunctions.hh"

#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Common/OSystem.hh"

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/StdTrsGeoBuilder.hh"

#include "CGNSWriter/CGNSWriter.hh"
#include "CGNSWriter/ParCGNSHighOrderWriter.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataHandleOutput.hh"

#include <unordered_set>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParCGNSHighOrderWriter, CGWriterData, CGNSWriterModule>
parCGNSHighOrderWriterProvider("ParCGNSHighOrderWriter");

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write CGNS file.");
}

//////////////////////////////////////////////////////////////////////////////

ParCGNSHighOrderWriter::ParCGNSHighOrderWriter(const std::string& name) : CGWriterCom(name),
    m_stdTrsGeoBuilder()
{
  addConfigOptionsTo(this);

  _fileFormatStr = "cgns";
  setParameter("FileFormat",&_fileFormatStr);

}
//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::execute()
{
  CFLog(INFO, "Writing solution to: " << getMethodData().getFilename().string() << "\n");

  writeToFile(getMethodData().getFilename().string());
  
}
//////////////////////////////////////////////////////////////////////////////

const std::string ParCGNSHighOrderWriter::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeToFile(const std::string& fileName)
{
  CFAUTOTRACE;

  // Initialization and MPI setup
  int rank, size;
  MPI_Comm comm = MPI_COMM_WORLD; // Use MPI_COMM_WORLD
  MPI_Info info = MPI_INFO_NULL;  // Use MPI_INFO_NULL

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Open the CGNS file in parallel
  openFile(fileName, comm, info);
  
  /*----------------------------! Retrieve general mesh and solution data -----------------------------------*/

  // get dimensionality, number of equations and reference length
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get the solution polynomial order
  const CFuint solOrder = (*elemType)[0].getSolOrder();
  // get the geometric order of the elements
  const CFuint geoOrder = (*elemType)[0].getGeoOrder(); 

  /*----------------------------! Define CGNS Base -----------------------------------*/

  const CFuint cellDim = dim; // Topological dimension of the cells (elements)
  const CFuint physDim = dim; // Physical dimension of the problem
  std::string BaseName = "Base1";

  // Create the CGNS base
  writeCGNSBase(cellDim, physDim, BaseName, fileIndex, baseIndex, comm);
  
  // Write time and iteration metadata
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFuint currentIter = subSysStatus->getNbIter();
  const CFreal currentTimeDim = subSysStatus->getCurrentTimeDim();
  const CFreal dt = subSysStatus->getDTDim();
  
  // Determine if this is a steady or unsteady case (dt <= 0 indicates steady)
  const bool isSteady = (dt <= 0.0);
  
  // Write simulation type (all ranks call this for consistency with parallel CGNS)
  CGNS_ENUMT(SimulationType_t) simType = isSteady ? CGNS_ENUMV(NonTimeAccurate) : CGNS_ENUMV(TimeAccurate);
  if (cg_simulation_type_write(fileIndex, baseIndex, simType) != CG_OK) {
    if (rank == 0) CFLog(VERBOSE, "CGNS Writer: Failed to write simulation type: " << cg_get_error() << "\n");
  }
  
  // Write BaseIterativeData (all ranks call this for consistency with parallel CGNS)
  int biter_result = cg_biter_write(fileIndex, baseIndex, "TimeIterValues", 1);
  if (biter_result == CG_OK) {
    cgsize_t nsteps = 1;
    int goto_result = cg_goto(fileIndex, baseIndex, "BaseIterativeData_t", 1, "end");
    if (goto_result == CG_OK) {
      // Write iteration number
      std::vector<int> iterValues{static_cast<int>(currentIter)};
      if (cg_array_write("IterationValues", Integer, 1, &nsteps, iterValues.data()) != CG_OK) {
        if (rank == 0) CFLog(VERBOSE, "CGNS Writer: Failed to write IterationValues: " << cg_get_error() << "\n");
      }
      
      // Write time values: iteration for steady, physical time for unsteady
      std::vector<double> timeValues{isSteady ? static_cast<double>(currentIter) : static_cast<double>(currentTimeDim)};
      if (cg_array_write("TimeValues", RealDouble, 1, &nsteps, timeValues.data()) != CG_OK) {
        if (rank == 0) CFLog(VERBOSE, "CGNS Writer: Failed to write TimeValues: " << cg_get_error() << "\n");
      }
    } else {
      if (rank == 0) CFLog(VERBOSE, "CGNS Writer: Failed to navigate to BaseIterativeData: " << cg_get_error() << "\n");
    }
  } else {
    if (rank == 0) CFLog(VERBOSE, "CGNS Writer: Failed to create BaseIterativeData: " << cg_get_error() << "\n");
  }

/*----------------------------! Get and clean nodes and connectivity -----------------------------------*/

  // 1) Retrieve raw nodes & connectivity
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes =
      MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  Common::SafePtr<Common::ConnectivityTable<CFuint>> cellNodesConn =
      MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // 2) Build a local connectivity subset that skips ghost elements
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  StdTrsGeoBuilder::GeoData& geoData = m_stdTrsGeoBuilder.getDataGE();
  geoData.trs = trs;

  // Data structures to store the "owned" connectivity
  std::vector<std::vector<CFuint>> localRawConn;

  CFuint startCellAccum = 0;

  // Build localRawConn in 0-based node indexing
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType) {
      // get the local # of cells for this type
      SafePtr<vector<ElementTypeData>> elemType =
          MeshDataStack::getActive()->getElementTypeData();
      const CFuint nbLocalCells = (*elemType)[iElemType].getNbElems();
      const CFuint startCellIdx = (*elemType)[iElemType].getStartIdx();

      for (CFuint iCellLocal = 0; iCellLocal < nbLocalCells; ++iCellLocal) {
          geoData.idx = startCellIdx + iCellLocal;
          GeometricEntity* cell = m_stdTrsGeoBuilder.buildGE();

          vector<State*>* cellStates = cell->getStates();
          bool isOwned = (!cellStates->empty() && (*cellStates)[0]->isParUpdatable());

          m_stdTrsGeoBuilder.releaseGE();

          if (!isOwned)
              continue; // skip ghost

          // gather connectivity from cellNodesConn
          const CFuint rowIdx = startCellAccum + iCellLocal;
          const CFuint nbNodesInCell = cellNodesConn->nbCols(rowIdx);

          std::vector<CFuint> thisCellConn;
          thisCellConn.reserve(nbNodesInCell);

          for (CFuint j = 0; j < nbNodesInCell; ++j) {
              // 0-based node index from the original connectivity
              CFuint nodeID = (*cellNodesConn)(rowIdx, j);
              thisCellConn.push_back(nodeID);
          }
          localRawConn.push_back(std::move(thisCellConn));
      }
      startCellAccum += nbLocalCells;
  }

  // 3) Create cleaned node data and connectivity
  std::vector<Framework::Node*> cleanedNodes;
  cleanedNodes.reserve(nodes.size()); // Reserve memory for efficiency

  for (CFuint i = 0; i < nodes.size(); ++i) {
      cleanedNodes.push_back(nodes[i]);
  }

  //initialise cleaned nodes to nodes
  std::vector<std::vector<CFuint>> cleanedConnectivity = localRawConn;
  // Build cleaned versions 
  cleanGhostData(cleanedNodes, cleanedConnectivity);

  // 5) Recompute local & global counts with the updated data
  const CFuint totNbNodes = cleanedNodes.size();
  const CFuint totNbCells = cleanedConnectivity.size();

  CFuint globalNbNodes = 0, globalNbCells = 0;
  MPI_Allreduce(&totNbNodes, &globalNbNodes, 1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&totNbCells, &globalNbCells, 1, MPI_UNSIGNED, MPI_SUM, comm);

  /*----------------------------! Mesh Geomtric Upgrade -----------------------------------*/

  //Geomtric order Q is capped at 4 (maximum supported by CGNS)
  const CFuint MaxGeoOrder = 4 ;
  const CFuint NewGeoOrder = std::max(std::min(MaxGeoOrder,solOrder), geoOrder);
  std::vector<std::vector<cgsize_t>> upgradedConnectivity(totNbCells); //totNbCells remains the same
  std::vector< RealVector > upgradedNodes;

  upgradeGeoOrder(cleanedNodes, cleanedConnectivity, NewGeoOrder, upgradedNodes, upgradedConnectivity);
  if (NewGeoOrder >= geoOrder) //if we upgrade mesh
  {
    cleanData(upgradedNodes, upgradedConnectivity);
  }

  CFuint NewtotNbNodes = upgradedNodes.size();
  CFuint nbrNodesPerElem = upgradedConnectivity[0].size();

  CFuint NewglobalNbNodes;

  // Use MPI_Allreduce to sum up the local numbers across all ranks
  MPI_Allreduce(&NewtotNbNodes, &NewglobalNbNodes, 1, MPI_UNSIGNED, MPI_SUM, comm);

  // Check for potential overflow at very high core counts
  if (NewglobalNbNodes > static_cast<CFuint>(std::numeric_limits<cgsize_t>::max())) {
    std::cerr << "Rank " << rank << " ERROR: Total node count exceeds CGNS limits. Consider using fewer cores or mesh partitioning." << std::endl;
    MPI_Abort(comm, 1);
  }

  /*----------------------------! Define CGNS Zone -----------------------------------*/
  // Create Zone
  cgsize_t totalNumNodes = NewglobalNbNodes; 
  cgsize_t totalNumCells = globalNbCells;
  cgsize_t Size[3] = {totalNumNodes, totalNumCells, 0};
  std::string ZoneName = "Zone1";

  writeCGNSZone(Size, ZoneName, fileIndex, baseIndex, index_zone);

  /*----------------------------! Define Mesh Data -----------------------------------*/

  /*! Define Mesh Data: Nodes*/

  CFuint localnodeCount = NewtotNbNodes;  // Number of nodes this rank will write
  CFuint startnode;

  MPI_Exscan(&localnodeCount, &startnode, 1, MPI_UNSIGNED, MPI_SUM, comm);
  if (rank == 0) {
      startnode = 1;  // First rank starts at index 1
  } else {
      startnode += 1;  // Adjust for 1-based indexing
  }  
  cgsize_t start_n = startnode; // Starting index of nodes (1-based for CGNS)
  cgsize_t end_n = startnode+localnodeCount-1; // Ending index, assuming consecutive numbering

  writeCGNSNodeCoordinates(NewtotNbNodes, dim, upgradedNodes, fileIndex, baseIndex, index_zone, index_coord_x, index_coord_y, index_coord_z, start_n, end_n);


  /*! Define Mesh Data: Connectivity*/

  // loop over elements Types 
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = totNbCells;//(*elemType)[iElemType].getNbElems();

    // get element shape
    const CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();    

    // Use getCGNSCellShape to determine the elementType for CGNS
    ElementType_t elementType = getCGNSCellShape(shape, NewGeoOrder);

    CFuint localElemCount = nbrElems;  // Number of elements this rank will write
    CFuint startElem;

    MPI_Exscan(&localElemCount, &startElem, 1, MPI_UNSIGNED, MPI_SUM, comm);
    if (rank == 0) {
        startElem = 1;  // First rank starts at index 1
    } else {
        startElem += 1;  // Adjust for 1-based indexing
    }  

    CFuint start_e = startElem; // Starting index of elements (1-based for CGNS)
    CFuint end_e = startElem+localElemCount-1; // Ending index, assuming consecutive numbering


    // create data node for elements
    if (cgp_section_write(fileIndex, baseIndex, index_zone,"ElemSection", elementType, 1, globalNbCells, 0, &index_section))
    {
      cgp_error_exit();
    }

    // Allocate memory for the connectivity array based on the number of elements and nodes per element
    cgsize_t* connectivity = new cgsize_t[nbrElems * nbrNodesPerElem];

    // Populate the connectivity array, adjusting indices according to CGNS conventions
    for (CFuint i = 0; i < nbrElems; ++i) {
        for (CFuint j = 0; j < nbrNodesPerElem; ++j) {
            // Assuming Connectivity is correctly ordered for CGNS, including 1-based indexing
            connectivity[i * nbrNodesPerElem + j] = upgradedConnectivity[i][j] + start_n - 1 ;
        }
    }

        /* write the element connectivity in parallel */
    if (cgp_elements_write_data(fileIndex, baseIndex, index_zone, index_section, start_e, end_e, connectivity))
    {
      cgp_error_exit();
    }

    // Clean up the dynamically allocated array
    delete[] connectivity;
  }

/*----------------------------! Define Solution Data -----------------------------------*/

  // get convective variable set
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  // get variable names
  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  //Below we get additional variables if any
  vector<std::string> extraVarNames;
  if (getMethodData().shouldPrintExtraValues()) {
    extraVarNames = updateVarSet->getExtraVarNames();
  }
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();
  std::vector< std::string > dh_varnames = datahandle_output->getVarNames();

  // Branch based on solution order: P0 uses CellCenter, P1+ uses Vertex
  const bool isP0 = (solOrder == 0);

  if (isP0) 
  {
    // P0: Write cell-centered solution
    writeCGNSSolutionNodeCellCenter("CellCenterSolution", fileIndex, baseIndex, index_zone, index_sol);
    
    // For P0, write cell-centered data directly (one value per cell)
    writeP0CellCenteredSolution(comm, rank, fileIndex, baseIndex, index_zone, index_sol,
                                 totNbCells, globalNbCells, nbEqs, varNames, extraVarNames,
                                 dh_varnames, updateVarSet, datahandle_output, trs, 
                                 geoData, elemType, nbrElemTypes);
  } 
  else 
  {
    // P1+: Write vertex solution 
    writeCGNSSolutionNode("VertexSolution", fileIndex, baseIndex, index_zone, index_sol);

    //Create Field data for this process 


    // prepares to loop over cells by getting the GeometricEntityPool
    // get inner cells TRS
    //SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
    //StdTrsGeoBuilder::GeoData& geoData = m_stdTrsGeoBuilder.getDataGE();
    geoData.trs = trs;

    // loop over elements Types 
    for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
    {
      // get element shape
      const CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();    

      // Use getCGNSCellShape to determine the elementType for CGNS
      ElementType_t elementType = getCGNSCellShape(shape, NewGeoOrder);

      // mapped coordinates of the solution points
      const vector< RealVector > outputPntsMappedCoords = getOutputPntsMappedCoords(elementType);

      // number of states in element type
      const CFuint nbrStates = (*elemType)[iElemType].getNbStates();

      // number of nodes in element type (original mesh)
      const CFuint nbrNodes = (*elemType)[iElemType].getNbNodes();

      // get the number of elements
      const CFuint nbrElems = totNbCells;//(*elemType)[iElemType].getNbElems();

      // get start index of this element type in global element list
      CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

      // evaluate the basis functions in the output points
      geoData.idx = cellIdx;
      GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();
      vector< RealVector > solShapeFuncs;
      // get the nodes
      vector<Node*>* cellNodes = cell->getNodes();
      cf_assert(cellNodes->size() == nbrNodes);

      CFuint nbrOutPnts = nbrNodesPerElem;
      for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
      {
        RealVector solShapeFunc =
            cell->computeShapeFunctionAtMappedCoord   (outputPntsMappedCoords[iPnt]);
        solShapeFuncs.push_back(solShapeFunc);
      }

      //release the GeometricEntity
      m_stdTrsGeoBuilder.releaseGE();

      // variable for solutions in output nodes
      vector< State > outputPntState;
      for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
      {
        RealVector aux2(nbEqs);
        State state(aux2,false);
        outputPntState.push_back(state);
      }

      // some helper states
      RealVector dimState(nbEqs);
      RealVector extraValues; // size will be set in the VarSet

      // Preparing arrays for each variable across all elements
      std::vector<std::vector<double>> variableDataArrays(nbEqs+ extraVarNames.size()+ dh_varnames.size(), std::vector<double>(NewglobalNbNodes));

      // Loop over elements
      for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
      {
        // build the GeometricEntity
        geoData.idx = cellIdx;
        GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();

        // get the nodes
        vector<Node*>* cellNodes = cell->getNodes();
        cf_assert(cellNodes->size() == nbrNodes);

        // get the states
        vector<State*>* cellStates = cell->getStates();
        cf_assert(cellStates->size() == nbrStates);


        // evaluate states at the output points
        for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
        {
          outputPntState[iPnt]  = 0.0;
          for (CFuint iState = 0; iState < nbrStates; ++iState)
          {  
            outputPntState[iPnt] += solShapeFuncs[iPnt][iState]*(*(*cellStates)[iState]);
          }
        }

        // Step 2: Populate variableDataArrays with solution data
        for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
        {
          CFuint GlobalID= upgradedConnectivity[iElem][iPnt] - 1 ;

          // get state in this point
          const RealVector& nodalState = outputPntState[iPnt];

          // Dimensionalize solution and compute extra values if necessary
          if (getMethodData().shouldPrintExtraValues()) {
              updateVarSet->setDimensionalValuesPlusExtraValues(nodalState, dimState, extraValues);
          } else {
              updateVarSet->setDimensionalValues(nodalState, dimState);
          }

          // Assigning state values to the corresponding position in the array
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              variableDataArrays[iEq][GlobalID] = dimState[iEq];
          }

          // Extra variables
          for (CFuint iExtra = 0; iExtra < extraVarNames.size(); ++iExtra) {
              variableDataArrays[nbEqs + iExtra][GlobalID] = extraValues[iExtra];
          }
              
          // datahandles with state based data
          for (CFuint iVar = 0; iVar < dh_varnames.size(); ++iVar)
          {
            DataHandleOutput::DataHandleInfo var_info = datahandle_output->getStateData(iVar);
            CFuint var_var = var_info.first;
            CFuint var_nbvars = var_info.second;
            DataHandle<CFreal> var = var_info.third;
            
            vector<CFreal> outputPntStateSockets;
            outputPntStateSockets.resize(nbrOutPnts);

            outputPntStateSockets[iPnt]  = 0.0;
            for (CFuint iState = 0; iState < nbrStates; ++iState)
            {
              outputPntStateSockets[iPnt] += solShapeFuncs[iPnt][iState]*var(((*cellStates)[iState])->getLocalID(), var_var, var_nbvars);
            }
            variableDataArrays[nbEqs + extraVarNames.size() + iVar][GlobalID] = outputPntStateSockets[iPnt];
          }
        }

        //release the GeometricEntity
        m_stdTrsGeoBuilder.releaseGE();
      }

      /* write the solution field data in parallel */
      cgsize_t start_s = start_n;
      cgsize_t end_s = end_n;

      // Add barrier to prevent I/O race conditions at very high core counts
      //MPI_Barrier(comm);

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) 
      {
        if (cgp_field_write(fileIndex, baseIndex, index_zone, index_sol, RealDouble, varNames[iEq].c_str(), &fieldIndex) != CG_OK) 
        {
          std::cerr << "Rank " << rank << " Error writing solution node for " << varNames[iEq] << ": " << cg_get_error() << std::endl;
        }      
        if (cgp_field_write_data(fileIndex, baseIndex, index_zone, index_sol, fieldIndex, &start_s, &end_s, variableDataArrays[iEq].data()) != CG_OK) {
          std::cerr << "Rank " << rank << " Error writing solution data for " << varNames[iEq] << ": " << cg_get_error() << std::endl;
        }
      }

      // Extra variables
      for (CFuint iExtra = 0; iExtra < extraVarNames.size(); ++iExtra) 
      {
        if (cgp_field_write(fileIndex, baseIndex, index_zone, index_sol, RealDouble, extraVarNames[iExtra].c_str(), &fieldIndex) != CG_OK) 
        {
          std::cerr << "Rank " << rank << " Error writing solution node for " << extraVarNames[iExtra] << ": " << cg_get_error() << std::endl;
        }

        if (cgp_field_write_data(fileIndex, baseIndex, index_zone, index_sol, fieldIndex, &start_s, &end_s, variableDataArrays[nbEqs + iExtra].data()) != CG_OK) {
          std::cerr << "Rank " << rank << " Error writing solution data for " << extraVarNames[iExtra] << ": " << cg_get_error() << std::endl;
        }
      }

      // Datahandles variables
      for (CFuint iDh = 0; iDh < dh_varnames.size(); ++iDh) 
      {
        if (cgp_field_write(fileIndex, baseIndex, index_zone, index_sol, RealDouble, dh_varnames[iDh].c_str(), &fieldIndex) != CG_OK) 
        {
          std::cerr << "Rank " << rank << " Error writing solution node for " << dh_varnames[iDh] << ": " << cg_get_error() << std::endl;
        }      
        if (cgp_field_write_data(fileIndex, baseIndex, index_zone, index_sol, fieldIndex, &start_s, &end_s, variableDataArrays[nbEqs + extraVarNames.size() + iDh].data()) != CG_OK) {
          std::cerr << "Rank " << rank << " Error writing solution data for " << dh_varnames[iDh] << ": " << cg_get_error() << std::endl;
        }
      }

    } // end P1+ vertex solution loop
  } // end if-else P0 vs P1+

  // Final barrier to ensure all I/O operations are complete before closing
  MPI_Barrier(comm);

  // Close the CGNS file
  if (cgp_close(fileIndex) != CG_OK) 
  {
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::cerr << "Rank " << rank << " Error closing CGNS file: " << cg_get_error() << std::endl;
    MPI_Abort(comm, 1);
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::openFile(const std::string& fileName, MPI_Comm comm, MPI_Info info)
{
  CFAUTOTRACE;

  if (cgp_open(fileName.c_str(), CG_MODE_WRITE, &fileIndex) != CG_OK) 
  {
    std::cerr << "Error opening CGNS file for writing: " << cg_get_error() << std::endl;
    MPI_Abort(comm, 1); // Abort if cannot open file
    return;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeCGNSBase(const CFuint cellDim, const CFuint physDim, std::string BaseName, int& fileIndex, int& baseIndex, MPI_Comm comm) 
{
  CFAUTOTRACE;
  // Create the CGNS base
  if (cg_base_write(fileIndex, BaseName.c_str(), cellDim, physDim, &baseIndex) != CG_OK) {
      std::cerr << "Error creating CGNS base node: " << cg_get_error() << std::endl;
      MPI_Abort(comm, 1); // Abort on severe error
      return; // Early return on failure
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeCGNSZone(cgsize_t size[3], std::string ZoneName, int& fileIndex, int& baseIndex, int& index_zone) 
{
  CFAUTOTRACE;

  if (cg_zone_write(fileIndex, baseIndex, ZoneName.c_str(), size, Unstructured, &index_zone) != CG_OK) {
      std::cerr << "Error creating unstructured zone: " << cg_get_error() << std::endl;
      cg_close(fileIndex); // Attempt to close the file on error
      return; // Early return on failure
  }

}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeCGNSNodeCoordinates(CFuint totNbNodes, CFuint dim, std::vector< RealVector > nodes,int& fileIndex, int& baseIndex, int& index_zone, int& index_coord_x, int& index_coord_y, int& index_coord_z,cgsize_t start,cgsize_t end) 
{
  CFAUTOTRACE;

  // create data nodes for coordinates
  if (cgp_coord_write(fileIndex, baseIndex, index_zone, RealDouble, "CoordinateX", &index_coord_x) ||
        cgp_coord_write(fileIndex, baseIndex, index_zone, RealDouble, "CoordinateY", &index_coord_y))
  {
    cgp_error_exit();
  }

  if (dim == 3) {
    if (cgp_coord_write(fileIndex, baseIndex, index_zone, RealDouble, "CoordinateZ", &index_coord_z))
    {
      cgp_error_exit();
    }
  }
        
  std::vector<double> xCoords(totNbNodes), yCoords(totNbNodes), zCoords(totNbNodes); // Use zCoords for 3D cases

  // Fill the coordinate arrays
  for (CFuint i = 0; i < totNbNodes; ++i) {
      xCoords[i] = nodes[i][0]; // Accessing x-coordinate
      yCoords[i] = nodes[i][1]; // Accessing y-coordinate
      if (dim == 3) {
          zCoords[i] = nodes[i][2]; // For 3D cases, also get the z-coordinate
      }
  }

  // Write coordinates to the CGNS file in parallel
  if (cgp_coord_write_data(fileIndex, baseIndex, index_zone, index_coord_x, &start, &end, xCoords.data()) ||
      cgp_coord_write_data(fileIndex, baseIndex, index_zone, index_coord_y, &start, &end, yCoords.data()))
  {
    cgp_error_exit();
  }
  if (dim == 3) {
    if (cgp_coord_write_data(fileIndex, baseIndex, index_zone, index_coord_z, &start, &end, zCoords.data()))
    {
      cgp_error_exit();
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeCGNSSolutionNode(const std::string& solutionName, int& fileIndex, int& baseIndex, int& index_zone, int& index_sol) 
{
  CFAUTOTRACE;

  // Attempt to write the solution node
  if ((cg_sol_write(fileIndex, baseIndex, index_zone, solutionName.c_str(), Vertex, &index_sol))!= CG_OK) 
  {
    std::cerr << "Error writing solution node: " << cg_get_error() << std::endl;
    cg_close(fileIndex); // Attempt to close the file on error
    return; 
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeCGNSSolutionNodeCellCenter(const std::string& solutionName, int& fileIndex, int& baseIndex, int& index_zone, int& index_sol) 
{
  CFAUTOTRACE;
  
  // Create cell-centered solution node (all ranks call this in parallel CGNS)
  if ((cg_sol_write(fileIndex, baseIndex, index_zone, solutionName.c_str(), CellCenter, &index_sol))!= CG_OK) 
  {
    std::cerr << "Error writing cell-centered solution node: " << cg_get_error() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::cleanData(std::vector<RealVector>& nodes, std::vector<std::vector<cgsize_t>>& connectivity) 
{
  CFAUTOTRACE;

  Common::CFMap<cgsize_t, cgsize_t> oldToNewIndexMap;
  std::vector<bool> isNodeUsed(nodes.size(), false);

  // Identify used nodes, considering the 1-based to 0-based indexing conversion
  for (const std::vector<cgsize_t>& cell : connectivity) {
      for (cgsize_t nodeId : cell) {
          // Convert from 1-based to 0-based indexing for internal representation
          cgsize_t zeroBasedNodeId = nodeId - 1;
          if (zeroBasedNodeId >= 0 && zeroBasedNodeId < isNodeUsed.size()) {
              isNodeUsed[zeroBasedNodeId] = true;
          } else {
              // Log warning or handle the case where nodeId is out of expected range
              CFLog(WARN, "Encountered nodeId " << nodeId << " out of range in connectivity. Please check your data.\n");
          }
      }
  }

  // Create a mapping for node indices
  cgsize_t newIndex = 0;
  for (cgsize_t i = 0; i < static_cast<cgsize_t>(isNodeUsed.size()); ++i) {
      if (isNodeUsed[i]) {
          oldToNewIndexMap.insert(i, newIndex++);
      }
  }

  oldToNewIndexMap.sortKeys(); // Sort the keys to prepare for lookups


  // Update connectivity with new indices, considering the 1-based to 0-based indexing conversion
  for (std::vector<cgsize_t>& cell : connectivity) {
      for (cgsize_t& nodeId : cell) {
          // Convert from 1-based to 0-based indexing
          cgsize_t zeroBasedNodeId = nodeId - 1;
          
          bool isFound;
          // Perform the lookup with the 0-based index
          cgsize_t newId = oldToNewIndexMap.find(zeroBasedNodeId, isFound);
          if (isFound) {
              // If found, convert back to 1-based indexing before storing
              nodeId = newId + 1;
          } else {
              // Handle the case where the key is not found, which might indicate an issue
              // or the fact that this nodeId was not used/updated
              CFLog(WARN, "Key not found in CFMap for nodeId " << nodeId << ". This nodeId may not be used in connectivity.\n");
          }
      }
  }

  // Remove unused nodes from the node list
  std::vector<RealVector> cleanedNodes;
  for (cgsize_t i = 0; i < static_cast<cgsize_t>(nodes.size()); ++i) {
      if (isNodeUsed[i]) {
          cleanedNodes.push_back(nodes[i]);
      }
  }
  nodes.swap(cleanedNodes); // Update the original nodes vector
}

//////////////////////////////////////////////////////////////////////////////
void ParCGNSHighOrderWriter::cleanGhostData(
    std::vector<Framework::Node*>& nodes,
    std::vector<std::vector<CFuint>>& connectivity)
{
    Common::CFMap<CFuint, CFuint> oldToNewIndexMap;
    std::vector<bool> isNodeUsed(nodes.size(), false);

    // Identify used nodes
    for (const auto& cell : connectivity) {
        for (CFuint nodeId : cell) {
            if (nodeId < nodes.size()) {
                isNodeUsed[nodeId] = true;
            } else {
                CFLog(INFO, "Encountered nodeId " << nodeId << " out of range in connectivity. Please check your data.\n");
            }
        }
    }

    // Create a mapping for node indices
    CFuint newIndex = 0;
    for (CFuint i = 0; i < static_cast<CFuint>(isNodeUsed.size()); ++i) {
        if (isNodeUsed[i]) {
            oldToNewIndexMap.insert(i, newIndex++);
        }
    }

    oldToNewIndexMap.sortKeys(); // Sort the keys for consistent lookups

    // Update connectivity with new indices
    for (auto& cell : connectivity) {
        for (CFuint& nodeId : cell) {
            bool isFound;
            CFuint newId = oldToNewIndexMap.find(nodeId, isFound);
            if (isFound) {
                nodeId = newId;
            } else {
                CFLog(INFO, "Key not found in CFMap for nodeId " << nodeId << ". This nodeId may not be used in connectivity.\n");
            }
        }
    }

    // Remove unused nodes from the node list
    std::vector<Framework::Node*> cleanedNodes;
    for (CFuint i = 0; i < static_cast<CFuint>(nodes.size()); ++i) {
        if (isNodeUsed[i]) {
            cleanedNodes.push_back(nodes[i]);
        }
    }
    nodes.swap(cleanedNodes); // Update the original nodes vector
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::setup()
{
  CFAUTOTRACE;

  CGWriterCom::setup();

  // setup geobuilder
  m_stdTrsGeoBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeBoundarySurface()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ParCGNSHighOrderWriter::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  return result;
}
//////////////////////////////////////////////////////////////////////////////

bool ParCGNSHighOrderWriter::nodesAreClose(const RealVector& a, const RealVector& b, int dim) {
    const double tol = 1e-8;
    for (int i = 0; i < dim; ++i) {
        if (std::abs(a[i] - b[i]) >= tol) return false;
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::upgradeGeoOrder(std::vector<Framework::Node*> nodes, std::vector<std::vector<CFuint>> cellNodesConn, CFuint NewGeoOrder, std::vector< RealVector >& upgradedNodes,  std::vector<std::vector<cgsize_t>>& upgradedConnectivity) 
{
  CFAUTOTRACE;

  // Initial setup: Retrieving existing mesh data
  const CFuint totNbNodes = nodes.size();
  const CFuint totNbCells = cellNodesConn.size();//cellNodesConn->nbRows();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");

  // Get the GeometricEntityPool
  StdTrsGeoBuilder::GeoData& geoData = m_stdTrsGeoBuilder.getDataGE();
  geoData.trs = trs;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // Step 1: Initialize upgradedNodes with existing node positions

  CFuint dim = (*(nodes)[0]).size();
  upgradedNodes.resize(totNbNodes);
  for (CFuint i = 0; i < totNbNodes; ++i) 
  {
    upgradedNodes[i].resize(dim);
    for (CFuint j = 0; j < dim; ++j) 
    {
      upgradedNodes[i][j]=((*(nodes)[i])[j]);
    }
  }

  // Step 2: Initialize the connectivity table with the old data but right CGNS ordering (if the initial mesh is Q2 take only corner nodes)
  CFuint geoOrder;
  // loop over elements Types 
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get element shape
    const CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();

    // get the geometric order of the element
    geoOrder = (*elemType)[iElemType].getGeoOrder(); 

    // Use getCGNSCellShape to determine the elementType for CGNS
    ElementType_t elementType = getCGNSCellShape(shape, geoOrder);

    // get the number of elements
    const CFuint nbrElems = totNbCells;//(*elemType)[iElemType].getNbElems();

    // number of nodes in element type
    const CFuint nbrNodes = (*elemType)[iElemType].getNbNodes();

    // get the number of corner nodes (Q1)
    const CFuint NbrNodesQ1 = (getOutputPntsMappedCoords(getCGNSCellShape(shape, 1))).size();

    // Determine the indices ordering following CGNS convention
    std::vector<CFuint> newIndices = getCGNSCellNodeConn(elementType);

    if (NewGeoOrder>geoOrder)
    {
      for (CFuint iElem = 0; iElem < nbrElems; ++iElem) 
      {
        upgradedConnectivity[iElem].resize(NbrNodesQ1);
        for (size_t j = 0; j < NbrNodesQ1; ++j) {
            // The original connectivity is directly usable, adjusted for 1-based indexing
            upgradedConnectivity[iElem][j] = (cellNodesConn)[iElem][newIndices[j]] + 1;//(*cellNodesConn)(iElem,newIndices[j]) + 1;
        }
      }
    }
    else
    {
      for (CFuint iElem = 0; iElem < nbrElems; ++iElem) 
      {
        upgradedConnectivity[iElem].resize(nbrNodes);
        for (size_t j = 0; j < nbrNodes; ++j) {
            // The original connectivity is directly usable, adjusted for 1-based indexing
            upgradedConnectivity[iElem][j] = (cellNodesConn)[iElem][newIndices[j]] + 1;//(*cellNodesConn)(iElem,newIndices[j]) + 1;
        }
      }
    }
  }

  //Upgrade phase
  if (NewGeoOrder>geoOrder)
  {
    // Step 3: Loop over elements again to upgrade their geometric order
    for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType) 
    {
      CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();
      // get the number of elements
      CFuint nbrElems = totNbCells;//(*elemType)[iElemType].getNbElems();
      // number of nodes in element type
      const CFuint nbrNodes = (*elemType)[iElemType].getNbNodes();
      // get start index of this element type in global element list
      CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

      // Use getCGNSCellShape to determine the new elementType for CGNS
      ElementType_t elementType = getCGNSCellShape(shape, NewGeoOrder);

      // get the number of corner nodes (Q1)
      const CFuint NbrNodesQ1 = (getOutputPntsMappedCoords(getCGNSCellShape(shape, 1))).size();

      // mapped coordinates of output points 
      const vector< RealVector >  outputPntsMappedCoords = getOutputPntsMappedCoords(elementType);

      // number of output points
      const CFuint nbrOutPnts = outputPntsMappedCoords.size();

      // evaluate the basis functions in the output points
      geoData.idx = cellIdx;
      GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();
      vector< RealVector > solShapeFuncs;
      // get the nodes
      vector<Node*>* cellNodes = cell->getNodes();
      cf_assert(cellNodes->size() == nbrNodes);

      vector< RealVector > geoShapeFuncs;
      for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
      {
        RealVector geoShapeFunc =
            cell->computeGeoShapeFunctionAtMappedCoord(outputPntsMappedCoords[iPnt]);
        geoShapeFuncs.push_back(geoShapeFunc);
      }
      //release the GeometricEntity
      m_stdTrsGeoBuilder.releaseGE();

      for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx) 
      {
        // build the GeometricEntity
        geoData.idx = cellIdx;
        GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();
        // get the nodes
        vector<Node*>* cellNodes = cell->getNodes();
        cf_assert(cellNodes->size() == nbrNodes);

        std::vector< RealVector > newNodePositions = calculateNewNodePositionsForElement(outputPntsMappedCoords, nbrNodes, cellNodes, geoShapeFuncs, NbrNodesQ1);

        CFuint startQ1Index = upgradedNodes.size();

        for (CFuint posIndex = 0; posIndex < newNodePositions.size(); ++posIndex) {
            RealVector newPos = newNodePositions[posIndex];
            bool nodeExists = false;
            CFuint nodeIndex = 0;

            // Start checking from the end of the initial Q1 nodes
            /*for (CFuint i = startQ1Index; i < upgradedNodes.size(); ++i) {
                if (nodesAreClose(upgradedNodes[i], newPos, dim)) {
                    nodeExists = true;
                    nodeIndex = i + 1; // Using 1-based indexing
                    break;
                }
            }

            if (!nodeExists) {
                upgradedNodes.push_back(newPos);
                nodeIndex = upgradedNodes.size(); // Using 1-based indexing
            }*/
                upgradedNodes.push_back(newPos);
                nodeIndex = upgradedNodes.size();
            // Add the index of the new node to upgradedConnectivity for this element
            upgradedConnectivity[iElem].push_back(nodeIndex);
        }

        //release the GeometricEntity
        m_stdTrsGeoBuilder.releaseGE();      
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

std::vector< RealVector > ParCGNSHighOrderWriter::calculateNewNodePositionsForElement(std::vector< RealVector >  outputPntsMappedCoords, CFuint nbrNodes, std::vector<COOLFluiD::Framework::Node*>*& cellNodes, std::vector< RealVector > geoShapeFuncs, CFuint NbrNodesQ1)
{
  CFAUTOTRACE;

  // number of output points
  const CFuint nbrOutPnts = outputPntsMappedCoords.size();
  // number of dims
  const CFuint dim = outputPntsMappedCoords[0].size();

  // variable for coordinates in output nodes
  std::vector< RealVector > outputPntCoords;
  outputPntCoords.resize(nbrOutPnts);
  for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
  {
    outputPntCoords[iPnt].resize(dim);
  }

  // evaluate node coordinates at the output points
  for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
  {
    outputPntCoords[iPnt] = 0.0;
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      for (CFuint iDim = 0; iDim < dim; ++iDim)
      {
        outputPntCoords[iPnt][iDim] += geoShapeFuncs[iPnt][iNode]*(*(*cellNodes)[iNode])[iDim];
      }
    }
  }

  // mapped coordinates of output points exluding corner nodes
  std::vector<RealVector>  result;
  for (size_t i = NbrNodesQ1; i < outputPntCoords.size(); ++i) {
            result.push_back(outputPntCoords[i]);
        }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

ElementType_t ParCGNSHighOrderWriter::getCGNSCellShape(CFGeoShape::Type shape,CFuint geoOrder)
{
  CFAUTOTRACE;

  switch (shape)
  {
    case CFGeoShape::TRIAG:
    {
      if (geoOrder==1)
      {
        return ElementType_t::TRI_3;
      }
      else if (geoOrder==2)
      {
        return ElementType_t::TRI_6;
      }
      else if (geoOrder==3)
      {
        return ElementType_t::TRI_10;
      }
      else if (geoOrder==4)
      {
        return ElementType_t::TRI_15;
      }
      break;
    }
    case CFGeoShape::QUAD:
    {
      if (geoOrder==1)
      {
        return ElementType_t::QUAD_4;
      }
      else if (geoOrder==2)
      {
        return ElementType_t::QUAD_9;
      }
      else if (geoOrder==3)
      {
        return ElementType_t::QUAD_16;
      }
      else if (geoOrder==4)
      {
        return ElementType_t::QUAD_25;
      }
      break;
    }
    case CFGeoShape::TETRA:
    {
      if (geoOrder==1)
      {
        return ElementType_t::TETRA_4;
      }
      else if (geoOrder==2)
      {
        return ElementType_t::TETRA_10;
      }
      else if (geoOrder==3)
      {
        return ElementType_t::TETRA_20;
      }
      else if (geoOrder==4)
      {
        return ElementType_t::TETRA_35;
      }
      break;
    }
    case CFGeoShape::PRISM:
    {
      if (geoOrder==1)
      {
        return ElementType_t::PENTA_6;
      }
      else if (geoOrder==2)
      {
        return ElementType_t::PENTA_18;
      }
      else if (geoOrder==3)
      {
        return ElementType_t::PENTA_40;
      }
      else if (geoOrder==4)
      {
        return ElementType_t::PENTA_75;
      }
      break;
    }
    case CFGeoShape::HEXA:
    {
      if (geoOrder==1)
      {
        return ElementType_t::HEXA_8;
      }
      else if (geoOrder==2)
      {
        return ElementType_t::HEXA_27;
      }
      else if (geoOrder==3)
      {
        return ElementType_t::HEXA_64;
      }
      else if (geoOrder==4)
      {
        return ElementType_t::HEXA_125;
      }
      break;
    }
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"ParCGNSHighOrderWriter::getCGNSCellShape Cell shape not supported or not implemented");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<CFuint> ParCGNSHighOrderWriter::getCGNSCellNodeConn(ElementType_t elementType)
{
  CFAUTOTRACE;

    switch (elementType)
    {
        case ElementType_t::BAR_2:
        {
          return std::vector<CFuint>{0, 1}; 
        }
        case ElementType_t::BAR_3:
        {
          return std::vector<CFuint>{0, 1, 2}; 
        }
        case ElementType_t::TRI_3:
        {
          return std::vector<CFuint>{0, 1, 2}; 
        }
        case ElementType_t::TRI_6:
        {
          return std::vector<CFuint>{0, 1, 2, 3, 4, 5}; 
        }
        case ElementType_t::QUAD_4:
        {
          return std::vector<CFuint>{0, 1, 2, 3};
        }
        case ElementType_t::QUAD_9:
        {
          return std::vector<CFuint>{0, 1, 2, 3, 4, 5, 6, 7, 8};
        }
        case ElementType_t::TETRA_4:
        {
          return std::vector<CFuint>{0, 1, 2, 3};
        }
        case ElementType_t::TETRA_10:
        {
            return std::vector<CFuint>{0, 1, 2, 3, 4, 5, 6, 9, 7, 8};
        }
        case ElementType_t::PENTA_6:
        {
          return std::vector<CFuint>{0, 1, 2, 3, 4, 5};
        }
        case ElementType_t::PENTA_18:
        {
            return std::vector<CFuint>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 16, 17, 10, 12, 14};
        }
        case ElementType_t::HEXA_8:
        {
          return std::vector<CFuint>{0, 1, 2, 3, 4, 5, 6, 7};
        }
        case ElementType_t::HEXA_27:
        {
            return std::vector<CFuint>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 22, 23, 24, 25, 12, 14, 16, 18, 20, 26, 21}; 
        }
        default:
        {
            throw Common::ShouldNotBeHereException(FromHere(),"Unsupported ElementType_t in getCGNSCellNodeConn");
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< RealVector > ParCGNSHighOrderWriter::getOutputPntsMappedCoords(ElementType_t type)
{
  CFAUTOTRACE;

  // mapped coordinates of output points 
  vector< RealVector >  outputPntsMappedCoords;

  switch (type)
  {
    case ElementType_t::TRI_3:
    {
      outputPntsMappedCoords.resize(3);
      for (CFuint i = 0; i < 3; ++i) {outputPntsMappedCoords[i].resize(2);}
      outputPntsMappedCoords[0][0] = 0.0;
      outputPntsMappedCoords[0][1] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; 
      outputPntsMappedCoords[1][1] = 0.0; 
      outputPntsMappedCoords[2][0] = 0.0; 
      outputPntsMappedCoords[2][1] = 1.0; 
      break;
    }
    case ElementType_t::TRI_6:
    {
      outputPntsMappedCoords.resize(6);
      for (CFuint i = 0; i < 6; ++i) {outputPntsMappedCoords[i].resize(2);}
      outputPntsMappedCoords[0][0] = 0.0;
      outputPntsMappedCoords[0][1] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; 
      outputPntsMappedCoords[1][1] = 0.0; 
      outputPntsMappedCoords[2][0] = 0.0; 
      outputPntsMappedCoords[2][1] = 1.0;  
      outputPntsMappedCoords[3][0] = 0.5; // Midpoint of edge 0-1
      outputPntsMappedCoords[3][1] = 0.0;
      outputPntsMappedCoords[4][0] = 0.5; // Midpoint of edge 1-2
      outputPntsMappedCoords[4][1] = 0.5;
      outputPntsMappedCoords[5][0] = 0.0; // Midpoint of edge 2-0
      outputPntsMappedCoords[5][1] = 0.5;
      break;
    }
    case ElementType_t::TRI_10:
    {
      outputPntsMappedCoords.resize(10);
      for (CFuint i = 0; i < 10; ++i) {outputPntsMappedCoords[i].resize(2);}
      outputPntsMappedCoords[0][0] = 0.0;
      outputPntsMappedCoords[0][1] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; 
      outputPntsMappedCoords[1][1] = 0.0; 
      outputPntsMappedCoords[2][0] = 0.0; 
      outputPntsMappedCoords[2][1] = 1.0; 
      outputPntsMappedCoords[3][0] = 1.0 / 3.0;
      outputPntsMappedCoords[3][1] = 0.0;
      outputPntsMappedCoords[4][0] = 2.0 / 3.0;
      outputPntsMappedCoords[4][1] = 0.0;
      outputPntsMappedCoords[5][0] = 2.0 / 3.0;
      outputPntsMappedCoords[5][1] = 1.0 / 3.0;
      outputPntsMappedCoords[6][0] = 1.0 / 3.0;
      outputPntsMappedCoords[6][1] = 2.0 / 3.0;
      outputPntsMappedCoords[7][0] = 0.0;
      outputPntsMappedCoords[7][1] = 2.0 / 3.0;
      outputPntsMappedCoords[8][0] = 0.0;
      outputPntsMappedCoords[8][1] = 1.0 / 3.0;
      outputPntsMappedCoords[9][0] = 1.0 / 3.0;
      outputPntsMappedCoords[9][1] = 1.0 / 3.0;
      break;
    }
    case ElementType_t::TRI_15:
    {
      outputPntsMappedCoords.resize(15);
      for (CFuint i = 0; i < 15; ++i) {outputPntsMappedCoords[i].resize(2);}
      outputPntsMappedCoords[0][0] = 0.0; // Node 1
      outputPntsMappedCoords[0][1] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; // Node 2
      outputPntsMappedCoords[1][1] = 0.0;
      outputPntsMappedCoords[2][0] = 0.0; // Node 3
      outputPntsMappedCoords[2][1] = 1.0;
      // Adding the edge nodes
      outputPntsMappedCoords[3][0] = 1.0 / 4.0; // Node 4
      outputPntsMappedCoords[3][1] = 0.0;
      outputPntsMappedCoords[4][0] = 2.0 / 4.0; // Node 5 (Midpoint of edge 1-2)
      outputPntsMappedCoords[4][1] = 0.0;
      outputPntsMappedCoords[5][0] = 3.0 / 4.0; // Node 6
      outputPntsMappedCoords[5][1] = 0.0;
      outputPntsMappedCoords[6][0] = 3.0 / 4.0; // Node 7
      outputPntsMappedCoords[6][1] = 1.0 / 4.0;
      outputPntsMappedCoords[7][0] = 2.0 / 4.0; // Node 8 (Midpoint of edge 2-3)
      outputPntsMappedCoords[7][1] = 2.0 / 4.0;
      outputPntsMappedCoords[8][0] = 1.0 / 4.0; // Node 9
      outputPntsMappedCoords[8][1] = 3.0 / 4.0;
      outputPntsMappedCoords[9][0] = 0.0; // Node 10
      outputPntsMappedCoords[9][1] = 3.0 / 4.0;
      outputPntsMappedCoords[10][0] = 0.0; // Node 11 (Midpoint of edge 3-1)
      outputPntsMappedCoords[10][1] = 2.0 / 4.0;
      outputPntsMappedCoords[11][0] = 0.0; // Node 12
      outputPntsMappedCoords[11][1] = 1.0 / 4.0;
      // Adding the interior nodes
      outputPntsMappedCoords[12][0] = 1.0 / 4.0; // Node 13
      outputPntsMappedCoords[12][1] = 1.0 / 4.0;
      outputPntsMappedCoords[13][0] = 2.0 / 4.0; // Node 14
      outputPntsMappedCoords[13][1] = 1.0 / 4.0;
      outputPntsMappedCoords[14][0] = 1.0 / 4.0; // Node 15
      outputPntsMappedCoords[14][1] = 2.0 / 4.0;
      break;
    }
    case ElementType_t::QUAD_4:
    {
      outputPntsMappedCoords.resize(4);
      for (CFuint i = 0; i < 4; ++i) {outputPntsMappedCoords[i].resize(2);}
      outputPntsMappedCoords[0][0]= -1.0;
      outputPntsMappedCoords[0][1]= -1.0;
      outputPntsMappedCoords[1][0]=  1.0; 
      outputPntsMappedCoords[1][1]= -1.0; 
      outputPntsMappedCoords[2][0]=  1.0; 
      outputPntsMappedCoords[2][1]=  1.0; 
      outputPntsMappedCoords[3][0]= -1.0; 
      outputPntsMappedCoords[3][1]=  1.0; 
      break;
    }
    case ElementType_t::QUAD_9:
    {
      outputPntsMappedCoords.resize(9);
      for (CFuint i = 0; i < 9; ++i) {outputPntsMappedCoords[i].resize(2);}
      outputPntsMappedCoords[0][0]= -1.0;
      outputPntsMappedCoords[0][1]= -1.0;
      outputPntsMappedCoords[1][0]=  1.0; 
      outputPntsMappedCoords[1][1]= -1.0; 
      outputPntsMappedCoords[2][0]=  1.0; 
      outputPntsMappedCoords[2][1]=  1.0; 
      outputPntsMappedCoords[3][0]= -1.0; 
      outputPntsMappedCoords[3][1]=  1.0; 

      outputPntsMappedCoords[4][0]=  0.0;
      outputPntsMappedCoords[4][1]= -1.0;
      outputPntsMappedCoords[5][0]=  1.0; 
      outputPntsMappedCoords[5][1]=  0.0; 
      outputPntsMappedCoords[6][0]=  0.0; 
      outputPntsMappedCoords[6][1]=  1.0; 
      outputPntsMappedCoords[7][0]= -1.0; 
      outputPntsMappedCoords[7][1]=  0.0;
      outputPntsMappedCoords[8][0]=  0.0; 
      outputPntsMappedCoords[8][1]=  0.0; 
      break;
    }
    case ElementType_t::QUAD_16:
    {
      outputPntsMappedCoords.resize(16);
      for (CFuint i = 0; i < 16; ++i) {outputPntsMappedCoords[i].resize(2);}
      // Corner nodes
      outputPntsMappedCoords[0][0] = -1.0; outputPntsMappedCoords[0][1] = -1.0;
      outputPntsMappedCoords[1][0] =  1.0; outputPntsMappedCoords[1][1] = -1.0;
      outputPntsMappedCoords[2][0] =  1.0; outputPntsMappedCoords[2][1] =  1.0;
      outputPntsMappedCoords[3][0] = -1.0; outputPntsMappedCoords[3][1] =  1.0;
      // Edge nodes
      // Bottom edge
      outputPntsMappedCoords[4][0] = -1.0/3.0; outputPntsMappedCoords[4][1] = -1.0;
      outputPntsMappedCoords[5][0] =  1.0/3.0; outputPntsMappedCoords[5][1] = -1.0;
      // Right edge
      outputPntsMappedCoords[6][0] = 1.0; outputPntsMappedCoords[6][1] = -1.0/3.0;
      outputPntsMappedCoords[7][0] = 1.0; outputPntsMappedCoords[7][1] =  1.0/3.0;
      // Top edge
      outputPntsMappedCoords[8][0] =  1.0/3.0; outputPntsMappedCoords[8][1] = 1.0;
      outputPntsMappedCoords[9][0] = -1.0/3.0; outputPntsMappedCoords[9][1] = 1.0;
      // Left edge
      outputPntsMappedCoords[10][0] = -1.0; outputPntsMappedCoords[10][1] =  1.0/3.0;
      outputPntsMappedCoords[11][0] = -1.0; outputPntsMappedCoords[11][1] = -1.0/3.0;
      // Face nodes
      outputPntsMappedCoords[12][0] = -1.0/3.0; outputPntsMappedCoords[12][1] = -1.0/3.0;
      outputPntsMappedCoords[13][0] =  1.0/3.0; outputPntsMappedCoords[13][1] = -1.0/3.0;
      outputPntsMappedCoords[14][0] =  1.0/3.0; outputPntsMappedCoords[14][1] =  1.0/3.0;
      outputPntsMappedCoords[15][0] = -1.0/3.0; outputPntsMappedCoords[15][1] =  1.0/3.0;
      break;
    }
    case ElementType_t::QUAD_25:
    {
      outputPntsMappedCoords.resize(25);
      for (CFuint i = 0; i < 25; ++i) {outputPntsMappedCoords[i].resize(2);}

      // Corner nodes
      outputPntsMappedCoords[0][0] = -1.0; outputPntsMappedCoords[0][1] = -1.0;
      outputPntsMappedCoords[1][0] =  1.0; outputPntsMappedCoords[1][1] = -1.0;
      outputPntsMappedCoords[2][0] =  1.0; outputPntsMappedCoords[2][1] =  1.0;
      outputPntsMappedCoords[3][0] = -1.0; outputPntsMappedCoords[3][1] =  1.0;

      // Edge nodes on the bottom edge
      outputPntsMappedCoords[4][0] = -0.5; outputPntsMappedCoords[4][1] = -1.0;
      outputPntsMappedCoords[5][0] = 0.0;  outputPntsMappedCoords[5][1] = -1.0;
      outputPntsMappedCoords[6][0] = 0.5;  outputPntsMappedCoords[6][1] = -1.0;
      // Right edge
      outputPntsMappedCoords[7][0] = 1.0; outputPntsMappedCoords[7][1] = -0.5;
      outputPntsMappedCoords[8][0] = 1.0;  outputPntsMappedCoords[8][1] = 0.0;
      outputPntsMappedCoords[9][0] = 1.0;  outputPntsMappedCoords[9][1] = 0.5;
      // Top edge
      outputPntsMappedCoords[10][0] = 0.5; outputPntsMappedCoords[10][1] = 1.0;
      outputPntsMappedCoords[11][0] = 0.0;  outputPntsMappedCoords[11][1] = 1.0;
      outputPntsMappedCoords[12][0] = -0.5;  outputPntsMappedCoords[12][1] = 1.0;
      // Left edge
      outputPntsMappedCoords[13][0] = -1.0; outputPntsMappedCoords[13][1] = 0.5;
      outputPntsMappedCoords[14][0] = -1.0;  outputPntsMappedCoords[14][1] = 0.0;
      outputPntsMappedCoords[15][0] = -1.0;  outputPntsMappedCoords[15][1] = -0.5;

      // Interior nodes in a 3x3 grid (row by row from bottom left corner)
      outputPntsMappedCoords[16][0] = -0.5; outputPntsMappedCoords[16][1] = -0.5;
      outputPntsMappedCoords[17][0] = 0.0;  outputPntsMappedCoords[17][1] = -0.5;
      outputPntsMappedCoords[18][0] = 0.5;  outputPntsMappedCoords[18][1] = -0.5;

      outputPntsMappedCoords[19][0] = 0.5; outputPntsMappedCoords[19][1] = 0.0;
      outputPntsMappedCoords[20][0] = 0.5;  outputPntsMappedCoords[20][1] = 0.5;
      outputPntsMappedCoords[21][0] = 0.0;  outputPntsMappedCoords[21][1] = 0.5;

      outputPntsMappedCoords[22][0] = -0.5; outputPntsMappedCoords[22][1] = 0.5;
      outputPntsMappedCoords[23][0] = -0.5;  outputPntsMappedCoords[23][1] = 0.0;
      outputPntsMappedCoords[24][0] = 0.0;  outputPntsMappedCoords[24][1] = 0.0;
      break;
    }
    case ElementType_t::TETRA_4:
    {
      outputPntsMappedCoords.resize(4);
      for (CFuint i = 0; i < 4; ++i) { outputPntsMappedCoords[i].resize(3); }
      // Vertices
      outputPntsMappedCoords[0][0] = 0.0; outputPntsMappedCoords[0][1] = 0.0; outputPntsMappedCoords[0][2] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; outputPntsMappedCoords[1][1] = 0.0; outputPntsMappedCoords[1][2] = 0.0;
      outputPntsMappedCoords[2][0] = 0.0; outputPntsMappedCoords[2][1] = 1.0; outputPntsMappedCoords[2][2] = 0.0;
      outputPntsMappedCoords[3][0] = 0.0; outputPntsMappedCoords[3][1] = 0.0; outputPntsMappedCoords[3][2] = 1.0;
      break;
    }
    case ElementType_t::TETRA_10:
    {
      outputPntsMappedCoords.resize(10);
      for (CFuint i = 0; i < 10; ++i) { outputPntsMappedCoords[i].resize(3); }

      // Vertices
      outputPntsMappedCoords[0][0] = 0.0; outputPntsMappedCoords[0][1] = 0.0; outputPntsMappedCoords[0][2] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; outputPntsMappedCoords[1][1] = 0.0; outputPntsMappedCoords[1][2] = 0.0;
      outputPntsMappedCoords[2][0] = 0.0; outputPntsMappedCoords[2][1] = 1.0; outputPntsMappedCoords[2][2] = 0.0;
      outputPntsMappedCoords[3][0] = 0.0; outputPntsMappedCoords[3][1] = 0.0; outputPntsMappedCoords[3][2] = 1.0;

      // Edge midpoints
      outputPntsMappedCoords[4][0] = 0.5; outputPntsMappedCoords[4][1] = 0.0; outputPntsMappedCoords[4][2] = 0.0; // Midpoint of edge 0-1
      outputPntsMappedCoords[5][0] = 0.5; outputPntsMappedCoords[5][1] = 0.5; outputPntsMappedCoords[5][2] = 0.0; // Midpoint of edge 1-2
      outputPntsMappedCoords[6][0] = 0.0; outputPntsMappedCoords[6][1] = 0.5; outputPntsMappedCoords[6][2] = 0.0; // Midpoint of edge 2-0
      outputPntsMappedCoords[7][0] = 0.0; outputPntsMappedCoords[7][1] = 0.0; outputPntsMappedCoords[7][2] = 0.5; // Midpoint of edge 0-3
      outputPntsMappedCoords[8][0] = 0.5; outputPntsMappedCoords[8][1] = 0.0; outputPntsMappedCoords[8][2] = 0.5; // Midpoint of edge 1-3
      outputPntsMappedCoords[9][0] = 0.0; outputPntsMappedCoords[9][1] = 0.5; outputPntsMappedCoords[9][2] = 0.5; // Midpoint of edge 2-3
      break;
    }
    case ElementType_t::TETRA_20:
    {
      outputPntsMappedCoords.resize(20);
      for (CFuint i = 0; i < 20; ++i) { outputPntsMappedCoords[i].resize(3); }
      // Vertices
      outputPntsMappedCoords[0][0] = 0.0; outputPntsMappedCoords[0][1] = 0.0; outputPntsMappedCoords[0][2] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; outputPntsMappedCoords[1][1] = 0.0; outputPntsMappedCoords[1][2] = 0.0;
      outputPntsMappedCoords[2][0] = 0.0; outputPntsMappedCoords[2][1] = 1.0; outputPntsMappedCoords[2][2] = 0.0;
      outputPntsMappedCoords[3][0] = 0.0; outputPntsMappedCoords[3][1] = 0.0; outputPntsMappedCoords[3][2] = 1.0;

      // Edge midpoints
      outputPntsMappedCoords[4][0] = 1.0 / 3.0 ; outputPntsMappedCoords[4][1] = 0.0; outputPntsMappedCoords[4][2] = 0.0; 
      outputPntsMappedCoords[5][0] = 2.0 / 3.0; outputPntsMappedCoords[5][1] = 0.0; outputPntsMappedCoords[5][2] = 0.0; 
      outputPntsMappedCoords[6][0] = 2.0 / 3.0; outputPntsMappedCoords[6][1] = 1.0 / 3.0; outputPntsMappedCoords[6][2] = 0.0; 
      outputPntsMappedCoords[7][0] = 1.0 / 3.0; outputPntsMappedCoords[7][1] = 2.0 / 3.0; outputPntsMappedCoords[7][2] = 0.0; 
      outputPntsMappedCoords[8][0] = 0.0; outputPntsMappedCoords[8][1] = 2.0 / 3.0; outputPntsMappedCoords[8][2] = 0.0; 
      outputPntsMappedCoords[9][0] = 0.0; outputPntsMappedCoords[9][1] = 1.0 / 3.0; outputPntsMappedCoords[9][2] = 0.0; 

      outputPntsMappedCoords[10][0] = 0.0; outputPntsMappedCoords[10][1] = 0.0; outputPntsMappedCoords[10][2] = 1.0 / 3.0; 
      outputPntsMappedCoords[11][0] = 0.0; outputPntsMappedCoords[11][1] = 0.0; outputPntsMappedCoords[11][2] = 2.0 / 3.0; 
      outputPntsMappedCoords[12][0] = 2.0 / 3.0; outputPntsMappedCoords[12][1] = 0.0; outputPntsMappedCoords[12][2] = 1.0 / 3.0; 
      outputPntsMappedCoords[13][0] = 1.0 / 3.0; outputPntsMappedCoords[13][1] = 0.0; outputPntsMappedCoords[13][2] = 2.0 / 3.0; 
      outputPntsMappedCoords[14][0] = 0.0; outputPntsMappedCoords[14][1] = 2.0 / 3.0; outputPntsMappedCoords[14][2] = 1.0 / 3.0; 
      outputPntsMappedCoords[15][0] = 0.0; outputPntsMappedCoords[15][1] = 1.0 / 3.0; outputPntsMappedCoords[15][2] = 2.0 / 3.0; 

      outputPntsMappedCoords[16][0] = 1.0 / 3.0; outputPntsMappedCoords[16][1] = 1.0 / 3.0; outputPntsMappedCoords[16][2] = 0.0; 
      outputPntsMappedCoords[17][0] = 1.0 / 3.0; outputPntsMappedCoords[17][1] = 0.0; outputPntsMappedCoords[17][2] = 1.0 / 3.0; 
      outputPntsMappedCoords[18][0] = 1.0 / 3.0; outputPntsMappedCoords[18][1] = 1.0 / 3.0; outputPntsMappedCoords[18][2] = 1.0 / 3.0; 
      outputPntsMappedCoords[19][0] = 0.0; outputPntsMappedCoords[19][1] = 1.0 / 3.0; outputPntsMappedCoords[19][2] = 1.0 / 3.0; 
      break;
    }
    case ElementType_t::TETRA_35:
    {
      outputPntsMappedCoords.resize(35);
      for (CFuint i = 0; i < 35; ++i) { outputPntsMappedCoords[i].resize(3); }
      
      // Vertices
      outputPntsMappedCoords[0][0] = 0.0; outputPntsMappedCoords[0][1] = 0.0; outputPntsMappedCoords[0][2] = 0.0;
      outputPntsMappedCoords[1][0] = 1.0; outputPntsMappedCoords[1][1] = 0.0; outputPntsMappedCoords[1][2] = 0.0;
      outputPntsMappedCoords[2][0] = 0.0; outputPntsMappedCoords[2][1] = 1.0; outputPntsMappedCoords[2][2] = 0.0;
      outputPntsMappedCoords[3][0] = 0.0; outputPntsMappedCoords[3][1] = 0.0; outputPntsMappedCoords[3][2] = 1.0;

      // Edge nodes at 0.25, 0.5, and 0.75 from each vertex
      // Edge 0-1
      outputPntsMappedCoords[4][0] = 0.25; outputPntsMappedCoords[4][1] = 0.0;  outputPntsMappedCoords[4][2] = 0.0;
      outputPntsMappedCoords[5][0] = 0.50; outputPntsMappedCoords[5][1] = 0.0;  outputPntsMappedCoords[5][2] = 0.0;
      outputPntsMappedCoords[6][0] = 0.75; outputPntsMappedCoords[6][1] = 0.0;  outputPntsMappedCoords[6][2] = 0.0;
      // Edge 1-2
      outputPntsMappedCoords[7][0] = 0.75; outputPntsMappedCoords[7][1] = 0.25;  outputPntsMappedCoords[7][2] = 0.0;
      outputPntsMappedCoords[8][0] = 0.50; outputPntsMappedCoords[8][1] = 0.5;  outputPntsMappedCoords[8][2] = 0.0;
      outputPntsMappedCoords[9][0] = 0.25; outputPntsMappedCoords[9][1] = 0.75;  outputPntsMappedCoords[9][2] = 0.0;
      // Edge 2-0
      outputPntsMappedCoords[10][0] = 0.0; outputPntsMappedCoords[10][1] = 0.75;  outputPntsMappedCoords[10][2] = 0.0;
      outputPntsMappedCoords[11][0] = 0.0; outputPntsMappedCoords[11][1] = 0.5;  outputPntsMappedCoords[11][2] = 0.0;
      outputPntsMappedCoords[12][0] = 0.0; outputPntsMappedCoords[12][1] = 0.25;  outputPntsMappedCoords[12][2] = 0.0;
      // Similar definitions for other edges
      outputPntsMappedCoords[13][0] = 0.0; outputPntsMappedCoords[13][1] = 0.0; outputPntsMappedCoords[13][2] = 0.25; 
      outputPntsMappedCoords[14][0] = 0.0; outputPntsMappedCoords[14][1] = 0.0; outputPntsMappedCoords[14][2] = 0.5 ; 
      outputPntsMappedCoords[15][0] = 0.0; outputPntsMappedCoords[15][1] = 0.0; outputPntsMappedCoords[15][2] = 0.75; 

      outputPntsMappedCoords[16][0] = 0.75; outputPntsMappedCoords[16][1] = 0.0; outputPntsMappedCoords[16][2] = 0.25; 
      outputPntsMappedCoords[17][0] = 0.5; outputPntsMappedCoords[17][1] = 0.0; outputPntsMappedCoords[17][2] = 0.5 ;
      outputPntsMappedCoords[18][0] = 0.25; outputPntsMappedCoords[18][1] = 0.0; outputPntsMappedCoords[18][2] = 0.75;
       
      outputPntsMappedCoords[19][0] = 0.0; outputPntsMappedCoords[19][1] = 0.75; outputPntsMappedCoords[19][2] = 0.25; 
      outputPntsMappedCoords[20][0] = 0.0; outputPntsMappedCoords[20][1] = 0.5 ; outputPntsMappedCoords[20][2] = 0.5 ; 
      outputPntsMappedCoords[21][0] = 0.0; outputPntsMappedCoords[21][1] = 0.25; outputPntsMappedCoords[21][2] = 0.75; 

      // Interior nodes 
      outputPntsMappedCoords[22][0] = 0.25; outputPntsMappedCoords[22][1] = 0.25; outputPntsMappedCoords[22][2] = 0.0; 
      outputPntsMappedCoords[23][0] = 0.5 ; outputPntsMappedCoords[23][1] = 0.25 ; outputPntsMappedCoords[23][2] = 0.0; 
      outputPntsMappedCoords[24][0] = 0.25; outputPntsMappedCoords[24][1] = 0.5; outputPntsMappedCoords[24][2] = 0.0; 

      // Face nodes - additional nodes on each face not on edges
      // These nodes form a triangle within each face of the tetrahedron
      outputPntsMappedCoords[25][0] = 0.25; outputPntsMappedCoords[25][1] = 0.0; outputPntsMappedCoords[25][2] = 0.25;
      outputPntsMappedCoords[26][0] = 0.50; outputPntsMappedCoords[26][1] = 0.0; outputPntsMappedCoords[26][2] = 0.25;

      outputPntsMappedCoords[27][0] = 0.25; outputPntsMappedCoords[27][1] = 0.0; outputPntsMappedCoords[27][2] = 0.50;

      outputPntsMappedCoords[28][0] = 0.50; outputPntsMappedCoords[28][1] = 0.25; outputPntsMappedCoords[28][2] = 0.25;
      outputPntsMappedCoords[29][0] = 0.25; outputPntsMappedCoords[29][1] = 0.5; outputPntsMappedCoords[29][2] = 0.25;

      outputPntsMappedCoords[30][0] = 0.25; outputPntsMappedCoords[30][1] = 0.25; outputPntsMappedCoords[30][2] = 0.5;

      outputPntsMappedCoords[31][0] = 0.; outputPntsMappedCoords[31][1] = 0.5; outputPntsMappedCoords[31][2] = 0.25;
      outputPntsMappedCoords[32][0] = 0.; outputPntsMappedCoords[32][1] = 0.25; outputPntsMappedCoords[32][2] = 0.25;

      outputPntsMappedCoords[33][0] = 0.; outputPntsMappedCoords[33][1] = 0.25; outputPntsMappedCoords[33][2] = 0.5;

      // Interior node - only one for TETRA_35, at the centroid
      outputPntsMappedCoords[34][0] = 0.25; outputPntsMappedCoords[34][1] = 0.25; outputPntsMappedCoords[34][2] = 0.25;
      break;
    }
    case ElementType_t::PENTA_6:
    {
      outputPntsMappedCoords.resize(6);
      for (CFuint i = 0; i < 6; ++i) {outputPntsMappedCoords[i].resize(3);}
      outputPntsMappedCoords[0][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[0][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[0][2] = -1.0; // z-coordinate

      outputPntsMappedCoords[1][0] = 1.0; // x-coordinate
      outputPntsMappedCoords[1][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[1][2] = -1.0; // z-coordinate

      outputPntsMappedCoords[2][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[2][1] = 1.0; // y-coordinate
      outputPntsMappedCoords[2][2] = -1.0; // z-coordinate

      outputPntsMappedCoords[3][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[3][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[3][2] = 1.0; // z-coordinate

      outputPntsMappedCoords[4][0] = 1.0; // x-coordinate
      outputPntsMappedCoords[4][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[4][2] = 1.0; // z-coordinate

      outputPntsMappedCoords[5][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[5][1] = 1.0; // y-coordinate
      outputPntsMappedCoords[5][2] = 1.0; // z-coordinate
      break;
    }
    case ElementType_t::PENTA_18:
    {
      outputPntsMappedCoords.resize(18);
      for (CFuint i = 0; i < 18; ++i) {outputPntsMappedCoords[i].resize(3);}
      outputPntsMappedCoords[0][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[0][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[0][2] = -1.0; // z-coordinate

      outputPntsMappedCoords[1][0] = 1.0; // x-coordinate
      outputPntsMappedCoords[1][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[1][2] = -1.0; // z-coordinate

      outputPntsMappedCoords[2][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[2][1] = 1.0; // y-coordinate
      outputPntsMappedCoords[2][2] = -1.0; // z-coordinate

      outputPntsMappedCoords[3][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[3][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[3][2] = 1.0; // z-coordinate

      outputPntsMappedCoords[4][0] = 1.0; // x-coordinate
      outputPntsMappedCoords[4][1] = 0.0; // y-coordinate
      outputPntsMappedCoords[4][2] = 1.0; // z-coordinate

      outputPntsMappedCoords[5][0] = 0.0; // x-coordinate
      outputPntsMappedCoords[5][1] = 1.0; // y-coordinate
      outputPntsMappedCoords[5][2] = 1.0; // z-coordinate

      // Edge midpoints on the triangular bases
      outputPntsMappedCoords[6][0] = 0.5; // Node 7
      outputPntsMappedCoords[6][1] = 0.0;
      outputPntsMappedCoords[6][2] = -1.0;

      outputPntsMappedCoords[7][0] = 0.5; // Node 8
      outputPntsMappedCoords[7][1] = 0.5;
      outputPntsMappedCoords[7][2] = -1.0;

      outputPntsMappedCoords[8][0] = 0.0; // Node 9
      outputPntsMappedCoords[8][1] = 0.5;
      outputPntsMappedCoords[8][2] = -1.0;

      // Edge midpoints on the quadrilateral sides and center
      outputPntsMappedCoords[9][0] = 0.0; // Node 10
      outputPntsMappedCoords[9][1] = 0.0;
      outputPntsMappedCoords[9][2] = 0.0;

      outputPntsMappedCoords[10][0] = 1.0; // Node 11
      outputPntsMappedCoords[10][1] = 0.0;
      outputPntsMappedCoords[10][2] = 0.0;

      outputPntsMappedCoords[11][0] = 0.0; // Node 12
      outputPntsMappedCoords[11][1] = 1.0;
      outputPntsMappedCoords[11][2] = 0.0;

      outputPntsMappedCoords[12][0] = 0.5; // Node 13
      outputPntsMappedCoords[12][1] = 0.0;
      outputPntsMappedCoords[12][2] = 1.0;

      outputPntsMappedCoords[13][0] = 0.5; // Node 14
      outputPntsMappedCoords[13][1] = 0.5;
      outputPntsMappedCoords[13][2] = 1.0;

      outputPntsMappedCoords[14][0] = 0.0; // Node 15
      outputPntsMappedCoords[14][1] = 0.5;
      outputPntsMappedCoords[14][2] = 1.0;

      // Additional edge midpoints
      outputPntsMappedCoords[15][0] = 0.5; // Node 16
      outputPntsMappedCoords[15][1] = 0.0;
      outputPntsMappedCoords[15][2] = 0.0;

      outputPntsMappedCoords[16][0] = 0.5; // Node 17
      outputPntsMappedCoords[16][1] = 0.5;
      outputPntsMappedCoords[16][2] = 0.0;

      outputPntsMappedCoords[17][0] = 0.0; // Node 18
      outputPntsMappedCoords[17][1] = 0.5;
      outputPntsMappedCoords[17][2] = 0.0;
      break;
    }
    case ElementType_t::PENTA_40:
    {
      outputPntsMappedCoords.resize(40);
      for (CFuint i = 0; i < 40; ++i) {outputPntsMappedCoords[i].resize(3);}

      // Base prism vertices, already defined in PENTA_6
      outputPntsMappedCoords[0][0] = 0.0; outputPntsMappedCoords[0][1] = 0.0; outputPntsMappedCoords[0][2] = -1.0;
      outputPntsMappedCoords[1][0] = 1.0; outputPntsMappedCoords[1][1] = 0.0; outputPntsMappedCoords[1][2] = -1.0;
      outputPntsMappedCoords[2][0] = 0.0; outputPntsMappedCoords[2][1] = 1.0; outputPntsMappedCoords[2][2] = -1.0;
      outputPntsMappedCoords[3][0] = 0.0; outputPntsMappedCoords[3][1] = 0.0; outputPntsMappedCoords[3][2] =  1.0;
      outputPntsMappedCoords[4][0] = 1.0; outputPntsMappedCoords[4][1] = 0.0; outputPntsMappedCoords[4][2] =  1.0;
      outputPntsMappedCoords[5][0] = 0.0; outputPntsMappedCoords[5][1] = 1.0; outputPntsMappedCoords[5][2] =  1.0;

      // Edges along the z = -1 face (base of the prism)
      outputPntsMappedCoords[6][0] = 1.0/3.0; outputPntsMappedCoords[6][1] = 0.0;    outputPntsMappedCoords[6][2] = -1.0;
      outputPntsMappedCoords[7][0] = 2.0/3.0; outputPntsMappedCoords[7][1] = 0.0;    outputPntsMappedCoords[7][2] = -1.0;
      outputPntsMappedCoords[8][0] = 2.0/3.0; outputPntsMappedCoords[8][1] = 1.0/3.0;    outputPntsMappedCoords[8][2] = -1.0;
      outputPntsMappedCoords[9][0] = 1.0/3.0; outputPntsMappedCoords[9][1] = 2.0/3.0;    outputPntsMappedCoords[9][2] = -1.0;
      outputPntsMappedCoords[10][0] = 0.0; outputPntsMappedCoords[10][1] = 2.0/3.0;    outputPntsMappedCoords[10][2] = -1.0;
      outputPntsMappedCoords[11][0] = 0.0; outputPntsMappedCoords[11][1] = 1.0/3.0;    outputPntsMappedCoords[11][2] = -1.0;
      
      //Vertical mid-edge vertices
      outputPntsMappedCoords[12][0] = 0.0; outputPntsMappedCoords[12][1] = 0.0; outputPntsMappedCoords[12][2] = -1.0/3.0;
      outputPntsMappedCoords[13][0] = 0.0; outputPntsMappedCoords[13][1] = 0.0; outputPntsMappedCoords[13][2] =  1.0/3.0;

      outputPntsMappedCoords[14][0] = 1.0; outputPntsMappedCoords[14][1] = 0.0; outputPntsMappedCoords[14][2] = -1.0/3.0;
      outputPntsMappedCoords[15][0] = 1.0; outputPntsMappedCoords[15][1] = 0.0; outputPntsMappedCoords[15][2] =  1.0/3.0;

      outputPntsMappedCoords[16][0] = 0.0; outputPntsMappedCoords[16][1] = 1.0; outputPntsMappedCoords[16][2] = -1.0/3.0;
      outputPntsMappedCoords[17][0] = 0.0; outputPntsMappedCoords[17][1] = 1.0; outputPntsMappedCoords[17][2] =  1.0/3.0;

      // Edges along the z = 1 face (base of the prism)
      outputPntsMappedCoords[18][0] = 1.0/3.0; outputPntsMappedCoords[18][1] = 0.0;    outputPntsMappedCoords[18][2] = 1.0;
      outputPntsMappedCoords[19][0] = 2.0/3.0; outputPntsMappedCoords[19][1] = 0.0;    outputPntsMappedCoords[19][2] = 1.0;
      outputPntsMappedCoords[20][0] = 2.0/3.0; outputPntsMappedCoords[20][1] = 1.0/3.0;    outputPntsMappedCoords[20][2] = 1.0;
      outputPntsMappedCoords[21][0] = 1.0/3.0; outputPntsMappedCoords[21][1] = 2.0/3.0;    outputPntsMappedCoords[21][2] = 1.0;
      outputPntsMappedCoords[22][0] = 0.0; outputPntsMappedCoords[22][1] = 2.0/3.0;    outputPntsMappedCoords[22][2] = 1.0;
      outputPntsMappedCoords[23][0] = 0.0; outputPntsMappedCoords[23][1] = 1.0/3.0;    outputPntsMappedCoords[23][2] = 1.0;

      outputPntsMappedCoords[24][0] = 1.0/3.0; outputPntsMappedCoords[24][1] = 1.0/3.0;    outputPntsMappedCoords[24][2] = -1.0;
      // ... other edge nodes for the base
      outputPntsMappedCoords[25][0] = 1.0/3.0; outputPntsMappedCoords[25][1] = 0.0;    outputPntsMappedCoords[25][2] = -1.0/3.0;
      outputPntsMappedCoords[26][0] = 2.0/3.0; outputPntsMappedCoords[26][1] = 0.0;    outputPntsMappedCoords[26][2] = -1.0/3.0;
      outputPntsMappedCoords[27][0] = 2.0/3.0; outputPntsMappedCoords[27][1] = 0.0;    outputPntsMappedCoords[27][2] = 1.0/3.0;
      outputPntsMappedCoords[28][0] = 1.0/3.0; outputPntsMappedCoords[28][1] = 0.0;    outputPntsMappedCoords[28][2] = 1.0/3.0;

      outputPntsMappedCoords[29][0] = 2.0/3.0; outputPntsMappedCoords[29][1] = 1.0/3.0;    outputPntsMappedCoords[29][2] = -1.0/3.0;
      outputPntsMappedCoords[30][0] = 1.0/3.0; outputPntsMappedCoords[30][1] = 2.0/3.0;    outputPntsMappedCoords[30][2] = -1.0/3.0;
      outputPntsMappedCoords[31][0] = 1.0/3.0; outputPntsMappedCoords[31][1] = 2.0/3.0;    outputPntsMappedCoords[31][2] = 1.0/3.0;
      outputPntsMappedCoords[32][0] = 2.0/3.0; outputPntsMappedCoords[32][1] = 1.0/3.0;    outputPntsMappedCoords[32][2] = 1.0/3.0;

      outputPntsMappedCoords[33][0] = 0.0; outputPntsMappedCoords[33][1] = 2.0/3.0;    outputPntsMappedCoords[33][2] = -1.0/3.0;
      outputPntsMappedCoords[34][0] = 0.0; outputPntsMappedCoords[34][1] = 1.0/3.0;    outputPntsMappedCoords[34][2] = -1.0/3.0;
      outputPntsMappedCoords[35][0] = 0.0; outputPntsMappedCoords[35][1] = 1.0/3.0;    outputPntsMappedCoords[35][2] = 1.0/3.0;
      outputPntsMappedCoords[36][0] = 0.0; outputPntsMappedCoords[36][1] = 2.0/3.0;    outputPntsMappedCoords[36][2] = 1.0/3.0;

      outputPntsMappedCoords[37][0] = 1.0/3.0; outputPntsMappedCoords[37][1] = 1.0/3.0;    outputPntsMappedCoords[37][2] = 1.0;

      outputPntsMappedCoords[38][0] = 1.0/3.0; outputPntsMappedCoords[38][1] = 1.0/3.0;    outputPntsMappedCoords[38][2] = -1.0/3.0;

      outputPntsMappedCoords[39][0] = 1.0/3.0; outputPntsMappedCoords[39][1] = 1.0/3.0;    outputPntsMappedCoords[39][2] = 1.0/3.0;
      break;
    }
    case ElementType_t::PENTA_75:
    {
      outputPntsMappedCoords.resize(75);
      for (CFuint i = 0; i < 75; ++i) {outputPntsMappedCoords[i].resize(3);}

      // Base prism vertices, already defined in PENTA_6
      outputPntsMappedCoords[0][0] = 0.0; outputPntsMappedCoords[0][1] = 0.0; outputPntsMappedCoords[0][2] = -1.0;
      outputPntsMappedCoords[1][0] = 1.0; outputPntsMappedCoords[1][1] = 0.0; outputPntsMappedCoords[1][2] = -1.0;
      outputPntsMappedCoords[2][0] = 0.0; outputPntsMappedCoords[2][1] = 1.0; outputPntsMappedCoords[2][2] = -1.0;
      outputPntsMappedCoords[3][0] = 0.0; outputPntsMappedCoords[3][1] = 0.0; outputPntsMappedCoords[3][2] =  1.0;
      outputPntsMappedCoords[4][0] = 1.0; outputPntsMappedCoords[4][1] = 0.0; outputPntsMappedCoords[4][2] =  1.0;
      outputPntsMappedCoords[5][0] = 0.0; outputPntsMappedCoords[5][1] = 1.0; outputPntsMappedCoords[5][2] =  1.0;

      // Edges along the z = -1 face (base of the prism)
      outputPntsMappedCoords[6][0] = 0.25; outputPntsMappedCoords[6][1] = 0.0;    outputPntsMappedCoords[6][2] = -1.0;
      outputPntsMappedCoords[7][0] = 0.50; outputPntsMappedCoords[7][1] = 0.0;    outputPntsMappedCoords[7][2] = -1.0;
      outputPntsMappedCoords[8][0] = 0.75; outputPntsMappedCoords[8][1] = 0.0;    outputPntsMappedCoords[8][2] = -1.0;

      outputPntsMappedCoords[9][0] = 0.75; outputPntsMappedCoords[9][1] = 0.25;    outputPntsMappedCoords[9][2] = -1.0;
      outputPntsMappedCoords[10][0] = 0.5; outputPntsMappedCoords[10][1] = 0.5;    outputPntsMappedCoords[10][2] = -1.0;
      outputPntsMappedCoords[11][0] = 0.25; outputPntsMappedCoords[11][1] = 0.75;    outputPntsMappedCoords[11][2] = -1.0;

      outputPntsMappedCoords[12][0] = 0.0; outputPntsMappedCoords[12][1] = 0.75;    outputPntsMappedCoords[12][2] = -1.0;
      outputPntsMappedCoords[13][0] = 0.0; outputPntsMappedCoords[13][1] = 0.5;    outputPntsMappedCoords[13][2] = -1.0;
      outputPntsMappedCoords[14][0] = 0.0; outputPntsMappedCoords[14][1] = 0.25;    outputPntsMappedCoords[14][2] = -1.0;
      
      //Vertical mid-edge vertices
      outputPntsMappedCoords[15][0] = 0.0; outputPntsMappedCoords[15][1] = 0.0; outputPntsMappedCoords[15][2] = -0.5;
      outputPntsMappedCoords[16][0] = 0.0; outputPntsMappedCoords[16][1] = 0.0; outputPntsMappedCoords[16][2] =  0.0;
      outputPntsMappedCoords[17][0] = 0.0; outputPntsMappedCoords[17][1] = 0.0; outputPntsMappedCoords[17][2] =  0.5;

      outputPntsMappedCoords[18][0] = 1.0; outputPntsMappedCoords[18][1] = 0.0; outputPntsMappedCoords[18][2] = -0.5;
      outputPntsMappedCoords[19][0] = 1.0; outputPntsMappedCoords[19][1] = 0.0; outputPntsMappedCoords[19][2] =  0.0;
      outputPntsMappedCoords[20][0] = 1.0; outputPntsMappedCoords[20][1] = 0.0; outputPntsMappedCoords[20][2] =  0.5;

      outputPntsMappedCoords[21][0] = 0.0; outputPntsMappedCoords[21][1] = 1.0; outputPntsMappedCoords[21][2] = -0.5;
      outputPntsMappedCoords[22][0] = 0.0; outputPntsMappedCoords[22][1] = 1.0; outputPntsMappedCoords[22][2] =  0.0;
      outputPntsMappedCoords[23][0] = 0.0; outputPntsMappedCoords[23][1] = 1.0; outputPntsMappedCoords[23][2] =  0.5;


      // Edges along the z = 1 face (base of the prism)
      outputPntsMappedCoords[24][0] = 0.25; outputPntsMappedCoords[24][1] = 0.0;    outputPntsMappedCoords[24][2] = 1.0;
      outputPntsMappedCoords[25][0] = 0.50; outputPntsMappedCoords[25][1] = 0.0;    outputPntsMappedCoords[25][2] = 1.0;
      outputPntsMappedCoords[26][0] = 0.75; outputPntsMappedCoords[26][1] = 0.0;    outputPntsMappedCoords[26][2] = 1.0;

      outputPntsMappedCoords[27][0] = 0.75; outputPntsMappedCoords[27][1] = 0.25;    outputPntsMappedCoords[27][2] = 1.0;
      outputPntsMappedCoords[28][0] = 0.5; outputPntsMappedCoords[28][1] = 0.5;    outputPntsMappedCoords[28][2] = 1.0;
      outputPntsMappedCoords[29][0] = 0.25; outputPntsMappedCoords[29][1] = 0.75;    outputPntsMappedCoords[29][2] = 1.0;

      outputPntsMappedCoords[30][0] = 0.0; outputPntsMappedCoords[30][1] = 0.75;    outputPntsMappedCoords[30][2] = 1.0;
      outputPntsMappedCoords[31][0] = 0.0; outputPntsMappedCoords[31][1] = 0.5;    outputPntsMappedCoords[31][2] = 1.0;
      outputPntsMappedCoords[32][0] = 0.0; outputPntsMappedCoords[32][1] = 0.25;    outputPntsMappedCoords[32][2] = 1.0;

      // Interior vertices along the z = -1 face (base of the prism)
      outputPntsMappedCoords[33][0] = 0.25; outputPntsMappedCoords[33][1] = 0.25;    outputPntsMappedCoords[33][2] = -1.0;
      outputPntsMappedCoords[34][0] = 0.5; outputPntsMappedCoords[34][1] = 0.25;    outputPntsMappedCoords[34][2] = -1.0;
      outputPntsMappedCoords[35][0] = 0.25; outputPntsMappedCoords[35][1] = 0.5;    outputPntsMappedCoords[35][2] = -1.0;
      // ... face nodes 1
      outputPntsMappedCoords[36][0] = 0.25; outputPntsMappedCoords[36][1] = 0.0;    outputPntsMappedCoords[36][2] = -0.5;
      outputPntsMappedCoords[37][0] = 0.50; outputPntsMappedCoords[37][1] = 0.0;    outputPntsMappedCoords[37][2] = -0.5;
      outputPntsMappedCoords[38][0] = 0.75; outputPntsMappedCoords[38][1] = 0.0;    outputPntsMappedCoords[38][2] = -0.5;
      outputPntsMappedCoords[39][0] = 0.75; outputPntsMappedCoords[39][1] = 0.0;    outputPntsMappedCoords[39][2] = 0.0;
      outputPntsMappedCoords[40][0] = 0.75; outputPntsMappedCoords[40][1] = 0.0;    outputPntsMappedCoords[40][2] = 0.5;
      outputPntsMappedCoords[41][0] = 0.5; outputPntsMappedCoords[41][1] = 0.0;    outputPntsMappedCoords[41][2] = 0.5;
      outputPntsMappedCoords[42][0] = 0.25; outputPntsMappedCoords[42][1] = 0.0;    outputPntsMappedCoords[42][2] = 0.5;
      outputPntsMappedCoords[43][0] = 0.25; outputPntsMappedCoords[43][1] = 0.0;    outputPntsMappedCoords[43][2] = 0.0;
      outputPntsMappedCoords[44][0] = 0.5; outputPntsMappedCoords[44][1] = 0.0;    outputPntsMappedCoords[44][2] = 0.0;
      // ... face nodes 2
      outputPntsMappedCoords[45][0] = 0.75; outputPntsMappedCoords[45][1] = 0.25;    outputPntsMappedCoords[45][2] = -0.5;
      outputPntsMappedCoords[46][0] = 0.5; outputPntsMappedCoords[46][1] = 0.5;    outputPntsMappedCoords[46][2] = -0.5;
      outputPntsMappedCoords[47][0] = 0.25; outputPntsMappedCoords[47][1] = 0.75;    outputPntsMappedCoords[47][2] = -0.5;
      outputPntsMappedCoords[48][0] = 0.25; outputPntsMappedCoords[48][1] = 0.75;    outputPntsMappedCoords[48][2] = 0.0;
      outputPntsMappedCoords[49][0] = 0.25; outputPntsMappedCoords[49][1] = 0.75;    outputPntsMappedCoords[49][2] = 0.5;
      outputPntsMappedCoords[50][0] = 0.5; outputPntsMappedCoords[50][1] = 0.5;    outputPntsMappedCoords[50][2] = 0.5;
      outputPntsMappedCoords[51][0] = 0.75; outputPntsMappedCoords[51][1] = 0.25;    outputPntsMappedCoords[51][2] = 0.5;
      outputPntsMappedCoords[52][0] = 0.75; outputPntsMappedCoords[52][1] = 0.25;    outputPntsMappedCoords[52][2] = 0.0;
      outputPntsMappedCoords[53][0] = 0.5; outputPntsMappedCoords[53][1] = 0.5;    outputPntsMappedCoords[53][2] = 0.0;
      // ... face nodes 3
      outputPntsMappedCoords[54][0] = 0.0; outputPntsMappedCoords[54][1] = 0.75;    outputPntsMappedCoords[54][2] = -0.5;
      outputPntsMappedCoords[55][0] = 0.0; outputPntsMappedCoords[55][1] = 0.50;    outputPntsMappedCoords[55][2] = -0.5;
      outputPntsMappedCoords[56][0] = 0.0; outputPntsMappedCoords[56][1] = 0.25;    outputPntsMappedCoords[56][2] = -0.5;
      outputPntsMappedCoords[57][0] = 0.0; outputPntsMappedCoords[57][1] = 0.25;    outputPntsMappedCoords[57][2] = 0.0;
      outputPntsMappedCoords[58][0] = 0.0; outputPntsMappedCoords[58][1] = 0.25;    outputPntsMappedCoords[58][2] = 0.5;
      outputPntsMappedCoords[59][0] = 0.0; outputPntsMappedCoords[59][1] = 0.50;    outputPntsMappedCoords[59][2] = 0.5;
      outputPntsMappedCoords[60][0] = 0.0; outputPntsMappedCoords[60][1] = 0.75;    outputPntsMappedCoords[60][2] = 0.5;
      outputPntsMappedCoords[61][0] = 0.0; outputPntsMappedCoords[61][1] = 0.75;    outputPntsMappedCoords[61][2] = 0.0;
      outputPntsMappedCoords[62][0] = 0.0; outputPntsMappedCoords[62][1] = 0.50;    outputPntsMappedCoords[62][2] = 0.0;
      // ... interior nodes along the z = 1 face
      outputPntsMappedCoords[63][0] = 0.25; outputPntsMappedCoords[63][1] = 0.25;    outputPntsMappedCoords[63][2] = 1.0;
      outputPntsMappedCoords[64][0] = 0.5; outputPntsMappedCoords[64][1] = 0.25;    outputPntsMappedCoords[64][2] = 1.0;
      outputPntsMappedCoords[65][0] = 0.25; outputPntsMappedCoords[65][1] = 0.5;    outputPntsMappedCoords[65][2] = 1.0;
      // ... interior nodes along the z = -0.5 face
      outputPntsMappedCoords[66][0] = 0.25; outputPntsMappedCoords[66][1] = 0.25;    outputPntsMappedCoords[66][2] = -0.5;
      outputPntsMappedCoords[67][0] = 0.5; outputPntsMappedCoords[67][1] = 0.25;    outputPntsMappedCoords[67][2] = -0.5;
      outputPntsMappedCoords[68][0] = 0.25; outputPntsMappedCoords[68][1] = 0.5;    outputPntsMappedCoords[68][2] = -0.5;
      // ... interior nodes along the z = 0.0 face
      outputPntsMappedCoords[69][0] = 0.25; outputPntsMappedCoords[69][1] = 0.25;    outputPntsMappedCoords[69][2] = 0.0;
      outputPntsMappedCoords[70][0] = 0.5; outputPntsMappedCoords[70][1] = 0.25;    outputPntsMappedCoords[70][2] = 0.0;
      outputPntsMappedCoords[71][0] = 0.25; outputPntsMappedCoords[71][1] = 0.5;    outputPntsMappedCoords[71][2] = 0.0;
      // ... interior nodes along the z = 0.5 face
      outputPntsMappedCoords[72][0] = 0.25; outputPntsMappedCoords[72][1] = 0.25;    outputPntsMappedCoords[72][2] = 0.5;
      outputPntsMappedCoords[73][0] = 0.5; outputPntsMappedCoords[73][1] = 0.25;    outputPntsMappedCoords[73][2] = 0.5;
      outputPntsMappedCoords[74][0] = 0.25; outputPntsMappedCoords[74][1] = 0.5;    outputPntsMappedCoords[74][2] = 0.5;
      break;
    }
    case ElementType_t::HEXA_8:
    {
      outputPntsMappedCoords.resize(8);
      for (CFuint i = 0; i < 8; ++i) { outputPntsMappedCoords[i].resize(3); }

      // Vertices
      outputPntsMappedCoords[0][0] = -1.0; outputPntsMappedCoords[0][1] = -1.0; outputPntsMappedCoords[0][2] = -1.0;
      outputPntsMappedCoords[1][0] =  1.0; outputPntsMappedCoords[1][1] = -1.0; outputPntsMappedCoords[1][2] = -1.0;
      outputPntsMappedCoords[2][0] =  1.0; outputPntsMappedCoords[2][1] =  1.0; outputPntsMappedCoords[2][2] = -1.0;
      outputPntsMappedCoords[3][0] = -1.0; outputPntsMappedCoords[3][1] =  1.0; outputPntsMappedCoords[3][2] = -1.0;
      outputPntsMappedCoords[4][0] = -1.0; outputPntsMappedCoords[4][1] = -1.0; outputPntsMappedCoords[4][2] =  1.0;
      outputPntsMappedCoords[5][0] =  1.0; outputPntsMappedCoords[5][1] = -1.0; outputPntsMappedCoords[5][2] =  1.0;
      outputPntsMappedCoords[6][0] =  1.0; outputPntsMappedCoords[6][1] =  1.0; outputPntsMappedCoords[6][2] =  1.0;
      outputPntsMappedCoords[7][0] = -1.0; outputPntsMappedCoords[7][1] =  1.0; outputPntsMappedCoords[7][2] =  1.0;

      break;
    }
    case ElementType_t::HEXA_27:
    {
      outputPntsMappedCoords.resize(27);
      for (CFuint i = 0; i < 27; ++i) { outputPntsMappedCoords[i].resize(3); }
      
      // Corner nodes
      outputPntsMappedCoords[0][0] = -1.0; outputPntsMappedCoords[0][1] = -1.0; outputPntsMappedCoords[0][2] = -1.0;
      outputPntsMappedCoords[1][0] =  1.0; outputPntsMappedCoords[1][1] = -1.0; outputPntsMappedCoords[1][2] = -1.0;
      outputPntsMappedCoords[2][0] =  1.0; outputPntsMappedCoords[2][1] =  1.0; outputPntsMappedCoords[2][2] = -1.0;
      outputPntsMappedCoords[3][0] = -1.0; outputPntsMappedCoords[3][1] =  1.0; outputPntsMappedCoords[3][2] = -1.0;
      outputPntsMappedCoords[4][0] = -1.0; outputPntsMappedCoords[4][1] = -1.0; outputPntsMappedCoords[4][2] =  1.0;
      outputPntsMappedCoords[5][0] =  1.0; outputPntsMappedCoords[5][1] = -1.0; outputPntsMappedCoords[5][2] =  1.0;
      outputPntsMappedCoords[6][0] =  1.0; outputPntsMappedCoords[6][1] =  1.0; outputPntsMappedCoords[6][2] =  1.0;
      outputPntsMappedCoords[7][0] = -1.0; outputPntsMappedCoords[7][1] =  1.0; outputPntsMappedCoords[7][2] =  1.0;

      // Mid-edge nodes
      outputPntsMappedCoords[8][0] =  0.0; outputPntsMappedCoords[8][1] = -1.0; outputPntsMappedCoords[8][2] = -1.0;
      outputPntsMappedCoords[9][0] =  1.0; outputPntsMappedCoords[9][1] =  0.0; outputPntsMappedCoords[9][2] = -1.0;
      outputPntsMappedCoords[10][0] = 0.0; outputPntsMappedCoords[10][1] =  1.0; outputPntsMappedCoords[10][2] = -1.0;
      outputPntsMappedCoords[11][0] = -1.0; outputPntsMappedCoords[11][1] =  0.0; outputPntsMappedCoords[11][2] = -1.0;
      outputPntsMappedCoords[12][0] = -1.0; outputPntsMappedCoords[12][1] = -1.0; outputPntsMappedCoords[12][2] =  0.0;
      outputPntsMappedCoords[13][0] =  1.0; outputPntsMappedCoords[13][1] = -1.0; outputPntsMappedCoords[13][2] =  0.0;
      outputPntsMappedCoords[14][0] =  1.0; outputPntsMappedCoords[14][1] =  1.0; outputPntsMappedCoords[14][2] =  0.0;
      outputPntsMappedCoords[15][0] = -1.0; outputPntsMappedCoords[15][1] =  1.0; outputPntsMappedCoords[15][2] =  0.0;
      outputPntsMappedCoords[16][0] =  0.0; outputPntsMappedCoords[16][1] = -1.0; outputPntsMappedCoords[16][2] =  1.0;
      outputPntsMappedCoords[17][0] =  1.0; outputPntsMappedCoords[17][1] =   0.0; outputPntsMappedCoords[17][2] =  1.0;
      outputPntsMappedCoords[18][0] =  0.0; outputPntsMappedCoords[18][1] =  1.0; outputPntsMappedCoords[18][2] =  1.0;
      outputPntsMappedCoords[19][0] = -1.0; outputPntsMappedCoords[19][1] =  0.0; outputPntsMappedCoords[19][2] =  1.0;

      // Face center nodes
      outputPntsMappedCoords[20][0] =  0.0; outputPntsMappedCoords[20][1] =  0.0; outputPntsMappedCoords[20][2] = -1.0;
      outputPntsMappedCoords[21][0] =  0.0; outputPntsMappedCoords[21][1] = -1.0; outputPntsMappedCoords[21][2] =  0.0;
      outputPntsMappedCoords[22][0] =  1.0; outputPntsMappedCoords[22][1] =  0.0; outputPntsMappedCoords[22][2] =  0.0;
      outputPntsMappedCoords[23][0] =  0.0; outputPntsMappedCoords[23][1] =  1.0; outputPntsMappedCoords[23][2] =  0.0;
      outputPntsMappedCoords[24][0] = -1.0; outputPntsMappedCoords[24][1] =  0.0; outputPntsMappedCoords[24][2] =  0.0;
      outputPntsMappedCoords[25][0] =  0.0; outputPntsMappedCoords[25][1] =  0.0; outputPntsMappedCoords[25][2] =  1.0;

      // Center node
      outputPntsMappedCoords[26][0] =  0.0; outputPntsMappedCoords[26][1] =  0.0; outputPntsMappedCoords[26][2] =  0.0;

      break;
    }
    case ElementType_t::HEXA_64:
    {
      outputPntsMappedCoords.resize(64);
      for (CFuint i = 0; i < 64; ++i) { outputPntsMappedCoords[i].resize(3); }
      // Corner nodes
      outputPntsMappedCoords[0][0] = -1.0; outputPntsMappedCoords[0][1] = -1.0; outputPntsMappedCoords[0][2] = -1.0;
      outputPntsMappedCoords[1][0] =  1.0; outputPntsMappedCoords[1][1] = -1.0; outputPntsMappedCoords[1][2] = -1.0;
      outputPntsMappedCoords[2][0] =  1.0; outputPntsMappedCoords[2][1] =  1.0; outputPntsMappedCoords[2][2] = -1.0;
      outputPntsMappedCoords[3][0] = -1.0; outputPntsMappedCoords[3][1] =  1.0; outputPntsMappedCoords[3][2] = -1.0;
      outputPntsMappedCoords[4][0] = -1.0; outputPntsMappedCoords[4][1] = -1.0; outputPntsMappedCoords[4][2] =  1.0;
      outputPntsMappedCoords[5][0] =  1.0; outputPntsMappedCoords[5][1] = -1.0; outputPntsMappedCoords[5][2] =  1.0;
      outputPntsMappedCoords[6][0] =  1.0; outputPntsMappedCoords[6][1] =  1.0; outputPntsMappedCoords[6][2] =  1.0;
      outputPntsMappedCoords[7][0] = -1.0; outputPntsMappedCoords[7][1] =  1.0; outputPntsMappedCoords[7][2] =  1.0;
      // Mid-edge nodes
      outputPntsMappedCoords[8][0] =  -1.0/3.0; outputPntsMappedCoords[8][1] = -1.0; outputPntsMappedCoords[8][2] = -1.0;
      outputPntsMappedCoords[9][0] =   1.0/3.0; outputPntsMappedCoords[9][1] = -1.0; outputPntsMappedCoords[9][2] = -1.0;
      
      outputPntsMappedCoords[10][0] = 1.0; outputPntsMappedCoords[10][1] =  -1.0/3.0; outputPntsMappedCoords[10][2] = -1.0;
      outputPntsMappedCoords[11][0] = 1.0; outputPntsMappedCoords[11][1] =   1.0/3.0; outputPntsMappedCoords[11][2] = -1.0;

      outputPntsMappedCoords[12][0] =  1.0/3.0; outputPntsMappedCoords[12][1] = 1.0; outputPntsMappedCoords[12][2] =  -1.0;
      outputPntsMappedCoords[13][0] = -1.0/3.0; outputPntsMappedCoords[13][1] = 1.0; outputPntsMappedCoords[13][2] =  -1.0;

      outputPntsMappedCoords[14][0] = -1.0; outputPntsMappedCoords[14][1] =  1.0/3.0; outputPntsMappedCoords[14][2] =  -1.0;
      outputPntsMappedCoords[15][0] = -1.0; outputPntsMappedCoords[15][1] = -1.0/3.0; outputPntsMappedCoords[15][2] =  -1.0;
      //---

      outputPntsMappedCoords[16][0] = -1.0; outputPntsMappedCoords[16][1] = -1.0; outputPntsMappedCoords[16][2] = -1.0/3.0;
      outputPntsMappedCoords[17][0] = -1.0; outputPntsMappedCoords[17][1] = -1.0; outputPntsMappedCoords[17][2] =  1.0/3.0;
      //---
      outputPntsMappedCoords[18][0] = 1.0; outputPntsMappedCoords[18][1] =  -1.0; outputPntsMappedCoords[18][2] = -1.0/3.0;
      outputPntsMappedCoords[19][0] = 1.0; outputPntsMappedCoords[19][1] =  -1.0; outputPntsMappedCoords[19][2] =  1.0/3.0;
      //---
      outputPntsMappedCoords[20][0] = 1.0; outputPntsMappedCoords[20][1] = 1.0; outputPntsMappedCoords[20][2] = -1.0/3.0;
      outputPntsMappedCoords[21][0] = 1.0; outputPntsMappedCoords[21][1] = 1.0; outputPntsMappedCoords[21][2] =  1.0/3.0;
      //---
      outputPntsMappedCoords[22][0] = -1.0; outputPntsMappedCoords[22][1] = 1.0; outputPntsMappedCoords[22][2] = -1.0/3.0;
      outputPntsMappedCoords[23][0] = -1.0; outputPntsMappedCoords[23][1] = 1.0; outputPntsMappedCoords[23][2] =  1.0/3.0;
      //---

      outputPntsMappedCoords[24][0] =  -1.0/3.0; outputPntsMappedCoords[24][1] = -1.0; outputPntsMappedCoords[24][2] = 1.0;
      outputPntsMappedCoords[25][0] =   1.0/3.0; outputPntsMappedCoords[25][1] = -1.0; outputPntsMappedCoords[25][2] = 1.0;
      
      outputPntsMappedCoords[26][0] = 1.0; outputPntsMappedCoords[26][1] =  -1.0/3.0; outputPntsMappedCoords[26][2] = 1.0;
      outputPntsMappedCoords[27][0] = 1.0; outputPntsMappedCoords[27][1] =   1.0/3.0; outputPntsMappedCoords[27][2] = 1.0;

      outputPntsMappedCoords[28][0] =  1.0/3.0; outputPntsMappedCoords[28][1] = 1.0; outputPntsMappedCoords[28][2] = 1.0;
      outputPntsMappedCoords[29][0] = -1.0/3.0; outputPntsMappedCoords[29][1] = 1.0; outputPntsMappedCoords[29][2] = 1.0;

      outputPntsMappedCoords[30][0] = -1.0; outputPntsMappedCoords[30][1] =  1.0/3.0; outputPntsMappedCoords[30][2] = 1.0;
      outputPntsMappedCoords[31][0] = -1.0; outputPntsMappedCoords[31][1] = -1.0/3.0; outputPntsMappedCoords[31][2] = 1.0;
      //---
      
      outputPntsMappedCoords[32][0] = -1.0/3.0; outputPntsMappedCoords[32][1] = -1.0/3.0;    outputPntsMappedCoords[32][2] = -1.0;
      outputPntsMappedCoords[33][0] =  1.0/3.0; outputPntsMappedCoords[33][1] = -1.0/3.0;    outputPntsMappedCoords[33][2] = -1.0;
      outputPntsMappedCoords[34][0] =  1.0/3.0; outputPntsMappedCoords[34][1] =  1.0/3.0;    outputPntsMappedCoords[34][2] = -1.0;
      outputPntsMappedCoords[35][0] = -1.0/3.0; outputPntsMappedCoords[35][1] =  1.0/3.0;    outputPntsMappedCoords[35][2] = -1.0;
      //---

      // ... face nodes 1
      outputPntsMappedCoords[36][0] = -1.0/3.0; outputPntsMappedCoords[36][1] = -1.0;    outputPntsMappedCoords[36][2] = -1.0/3.0;
      outputPntsMappedCoords[37][0] =  1.0/3.0; outputPntsMappedCoords[37][1] = -1.0;    outputPntsMappedCoords[37][2] = -1.0/3.0;
      outputPntsMappedCoords[38][0] =  1.0/3.0; outputPntsMappedCoords[38][1] = -1.0;    outputPntsMappedCoords[38][2] =  1.0/3.0;
      outputPntsMappedCoords[39][0] = -1.0/3.0; outputPntsMappedCoords[39][1] = -1.0;    outputPntsMappedCoords[39][2] =  1.0/3.0;
      // ... face nodes 2
      outputPntsMappedCoords[40][0] = 1.0; outputPntsMappedCoords[40][1] = -1.0/3.0;    outputPntsMappedCoords[40][2] = -1.0/3.0;
      outputPntsMappedCoords[41][0] = 1.0; outputPntsMappedCoords[41][1] =  1.0/3.0;    outputPntsMappedCoords[41][2] = -1.0/3.0;
      outputPntsMappedCoords[42][0] = 1.0; outputPntsMappedCoords[42][1] =  1.0/3.0;    outputPntsMappedCoords[42][2] =  1.0/3.0;
      outputPntsMappedCoords[43][0] = 1.0; outputPntsMappedCoords[43][1] = -1.0/3.0;    outputPntsMappedCoords[43][2] =  1.0/3.0;
      // ... face nodes 3
      outputPntsMappedCoords[44][0] =  1.0/3.0; outputPntsMappedCoords[44][1] = 1.0;    outputPntsMappedCoords[44][2] = -1.0/3.0;
      outputPntsMappedCoords[45][0] = -1.0/3.0; outputPntsMappedCoords[45][1] = 1.0;    outputPntsMappedCoords[45][2] = -1.0/3.0;
      outputPntsMappedCoords[46][0] = -1.0/3.0; outputPntsMappedCoords[46][1] = 1.0;    outputPntsMappedCoords[46][2] =  1.0/3.0;
      outputPntsMappedCoords[47][0] =  1.0/3.0; outputPntsMappedCoords[47][1] = 1.0;    outputPntsMappedCoords[47][2] =  1.0/3.0;
      // ... face nodes 4
      outputPntsMappedCoords[48][0] = -1.0; outputPntsMappedCoords[48][1] =  1.0/3.0;    outputPntsMappedCoords[48][2] = -1.0/3.0;
      outputPntsMappedCoords[49][0] = -1.0; outputPntsMappedCoords[49][1] = -1.0/3.0;    outputPntsMappedCoords[49][2] = -1.0/3.0;
      outputPntsMappedCoords[50][0] = -1.0; outputPntsMappedCoords[50][1] = -1.0/3.0;    outputPntsMappedCoords[50][2] =  1.0/3.0;
      outputPntsMappedCoords[51][0] = -1.0; outputPntsMappedCoords[51][1] =  1.0/3.0;    outputPntsMappedCoords[51][2] =  1.0/3.0;
      //---
      
      outputPntsMappedCoords[52][0] = -1.0/3.0; outputPntsMappedCoords[52][1] = -1.0/3.0;    outputPntsMappedCoords[52][2] = 1.0;
      outputPntsMappedCoords[53][0] =  1.0/3.0; outputPntsMappedCoords[53][1] = -1.0/3.0;    outputPntsMappedCoords[53][2] = 1.0;
      outputPntsMappedCoords[54][0] =  1.0/3.0; outputPntsMappedCoords[54][1] =  1.0/3.0;    outputPntsMappedCoords[54][2] = 1.0;
      outputPntsMappedCoords[55][0] = -1.0/3.0; outputPntsMappedCoords[55][1] =  1.0/3.0;    outputPntsMappedCoords[55][2] = 1.0;
      //---      
      
      outputPntsMappedCoords[56][0] = -1.0/3.0; outputPntsMappedCoords[56][1] = -1.0/3.0;    outputPntsMappedCoords[56][2] = -1.0/3.0;
      outputPntsMappedCoords[57][0] =  1.0/3.0; outputPntsMappedCoords[57][1] = -1.0/3.0;    outputPntsMappedCoords[57][2] = -1.0/3.0;
      outputPntsMappedCoords[58][0] =  1.0/3.0; outputPntsMappedCoords[58][1] =  1.0/3.0;    outputPntsMappedCoords[58][2] = -1.0/3.0;
      outputPntsMappedCoords[59][0] = -1.0/3.0; outputPntsMappedCoords[59][1] =  1.0/3.0;    outputPntsMappedCoords[59][2] = -1.0/3.0;
      //---          
      
      outputPntsMappedCoords[60][0] = -1.0/3.0; outputPntsMappedCoords[60][1] = -1.0/3.0;    outputPntsMappedCoords[60][2] = 1.0/3.0;
      outputPntsMappedCoords[61][0] =  1.0/3.0; outputPntsMappedCoords[61][1] = -1.0/3.0;    outputPntsMappedCoords[61][2] = 1.0/3.0;
      outputPntsMappedCoords[62][0] =  1.0/3.0; outputPntsMappedCoords[62][1] =  1.0/3.0;    outputPntsMappedCoords[62][2] = 1.0/3.0;
      outputPntsMappedCoords[63][0] = -1.0/3.0; outputPntsMappedCoords[63][1] =  1.0/3.0;    outputPntsMappedCoords[63][2] = 1.0/3.0;
      break;
    }
    case ElementType_t::HEXA_125:
    {
      outputPntsMappedCoords.resize(125);
      for (CFuint i = 0; i < 125; ++i) { outputPntsMappedCoords[i].resize(3); }
      outputPntsMappedCoords[0][0] = -1.0; outputPntsMappedCoords[0][1] = -1.0; outputPntsMappedCoords[0][2] = -1.0;
      outputPntsMappedCoords[1][0] =  1.0; outputPntsMappedCoords[1][1] = -1.0; outputPntsMappedCoords[1][2] = -1.0;
      outputPntsMappedCoords[2][0] =  1.0; outputPntsMappedCoords[2][1] =  1.0; outputPntsMappedCoords[2][2] = -1.0;
      outputPntsMappedCoords[3][0] = -1.0; outputPntsMappedCoords[3][1] =  1.0; outputPntsMappedCoords[3][2] = -1.0;
      outputPntsMappedCoords[4][0] = -1.0; outputPntsMappedCoords[4][1] = -1.0; outputPntsMappedCoords[4][2] =  1.0;
      outputPntsMappedCoords[5][0] =  1.0; outputPntsMappedCoords[5][1] = -1.0; outputPntsMappedCoords[5][2] =  1.0;
      outputPntsMappedCoords[6][0] =  1.0; outputPntsMappedCoords[6][1] =  1.0; outputPntsMappedCoords[6][2] =  1.0;
      outputPntsMappedCoords[7][0] = -1.0; outputPntsMappedCoords[7][1] =  1.0; outputPntsMappedCoords[7][2] =  1.0;
      // Mid-edge nodes
      outputPntsMappedCoords[8][0] = -0.5; outputPntsMappedCoords[8][1] = -1.0; outputPntsMappedCoords[8][2] = -1.0;
      outputPntsMappedCoords[9][0] =  0.0; outputPntsMappedCoords[9][1] = -1.0; outputPntsMappedCoords[9][2] = -1.0;
      outputPntsMappedCoords[10][0] =  0.5; outputPntsMappedCoords[10][1] = -1.0; outputPntsMappedCoords[10][2] = -1.0;
      
      outputPntsMappedCoords[11][0] = 1.0; outputPntsMappedCoords[11][1] = -0.5; outputPntsMappedCoords[11][2] = -1.0;
      outputPntsMappedCoords[12][0] = 1.0; outputPntsMappedCoords[12][1] =  0.0; outputPntsMappedCoords[12][2] = -1.0;
      outputPntsMappedCoords[13][0] = 1.0; outputPntsMappedCoords[13][1] =  0.5; outputPntsMappedCoords[13][2] = -1.0;

      outputPntsMappedCoords[14][0] =  0.5; outputPntsMappedCoords[14][1] = 1.0; outputPntsMappedCoords[14][2] =  -1.0;
      outputPntsMappedCoords[15][0] =  0.0; outputPntsMappedCoords[15][1] = 1.0; outputPntsMappedCoords[15][2] =  -1.0;
      outputPntsMappedCoords[16][0] = -0.5; outputPntsMappedCoords[16][1] = 1.0; outputPntsMappedCoords[16][2] =  -1.0;

      outputPntsMappedCoords[17][0] = -1.0; outputPntsMappedCoords[17][1] =  0.5; outputPntsMappedCoords[17][2] =  -1.0;
      outputPntsMappedCoords[18][0] = -1.0; outputPntsMappedCoords[18][1] =  0.0; outputPntsMappedCoords[18][2] =  -1.0;
      outputPntsMappedCoords[19][0] = -1.0; outputPntsMappedCoords[19][1] = -0.5; outputPntsMappedCoords[19][2] =  -1.0;
      //---

      outputPntsMappedCoords[20][0] = -1.0; outputPntsMappedCoords[20][1] = -1.0; outputPntsMappedCoords[20][2] = -0.5;
      outputPntsMappedCoords[21][0] = -1.0; outputPntsMappedCoords[21][1] = -1.0; outputPntsMappedCoords[21][2] =  0.0;
      outputPntsMappedCoords[22][0] = -1.0; outputPntsMappedCoords[22][1] = -1.0; outputPntsMappedCoords[22][2] =  0.5;
      //---
      outputPntsMappedCoords[23][0] = 1.0; outputPntsMappedCoords[23][1] =  -1.0; outputPntsMappedCoords[23][2] = -0.5;
      outputPntsMappedCoords[24][0] = 1.0; outputPntsMappedCoords[24][1] =  -1.0; outputPntsMappedCoords[24][2] =  0.0;
      outputPntsMappedCoords[25][0] = 1.0; outputPntsMappedCoords[25][1] =  -1.0; outputPntsMappedCoords[25][2] =  0.5;
      //---
      outputPntsMappedCoords[26][0] = 1.0; outputPntsMappedCoords[26][1] = 1.0; outputPntsMappedCoords[26][2] = -0.5;
      outputPntsMappedCoords[27][0] = 1.0; outputPntsMappedCoords[27][1] = 1.0; outputPntsMappedCoords[27][2] =  0.0;
      outputPntsMappedCoords[28][0] = 1.0; outputPntsMappedCoords[28][1] = 1.0; outputPntsMappedCoords[28][2] =  0.5;
      //---
      outputPntsMappedCoords[29][0] = -1.0; outputPntsMappedCoords[29][1] = 1.0; outputPntsMappedCoords[29][2] = -0.5;
      outputPntsMappedCoords[30][0] = -1.0; outputPntsMappedCoords[30][1] = 1.0; outputPntsMappedCoords[30][2] =  0.0;
      outputPntsMappedCoords[31][0] = -1.0; outputPntsMappedCoords[31][1] = 1.0; outputPntsMappedCoords[31][2] =  0.5;
      //---

      outputPntsMappedCoords[32][0] = -0.5; outputPntsMappedCoords[32][1] = -1.0; outputPntsMappedCoords[32][2] = 1.0;
      outputPntsMappedCoords[33][0] =  0.0; outputPntsMappedCoords[33][1] = -1.0; outputPntsMappedCoords[33][2] = 1.0;
      outputPntsMappedCoords[34][0] =  0.5; outputPntsMappedCoords[34][1] = -1.0; outputPntsMappedCoords[34][2] = 1.0;
      
      outputPntsMappedCoords[35][0] = 1.0; outputPntsMappedCoords[35][1] = -0.5; outputPntsMappedCoords[35][2] = 1.0;
      outputPntsMappedCoords[36][0] = 1.0; outputPntsMappedCoords[36][1] =  0.0; outputPntsMappedCoords[36][2] = 1.0;
      outputPntsMappedCoords[37][0] = 1.0; outputPntsMappedCoords[37][1] =  0.5; outputPntsMappedCoords[37][2] = 1.0;

      outputPntsMappedCoords[38][0] =  0.5; outputPntsMappedCoords[38][1] = 1.0; outputPntsMappedCoords[38][2] = 1.0;
      outputPntsMappedCoords[39][0] =  0.0; outputPntsMappedCoords[39][1] = 1.0; outputPntsMappedCoords[39][2] = 1.0;
      outputPntsMappedCoords[40][0] = -0.5; outputPntsMappedCoords[40][1] = 1.0; outputPntsMappedCoords[40][2] = 1.0;

      outputPntsMappedCoords[41][0] = -1.0; outputPntsMappedCoords[41][1] =  0.5; outputPntsMappedCoords[41][2] = 1.0;
      outputPntsMappedCoords[42][0] = -1.0; outputPntsMappedCoords[42][1] =  0.0; outputPntsMappedCoords[42][2] = 1.0;
      outputPntsMappedCoords[43][0] = -1.0; outputPntsMappedCoords[43][1] = -0.5; outputPntsMappedCoords[43][2] = 1.0;
      //---
      
      outputPntsMappedCoords[44][0] = -0.5; outputPntsMappedCoords[44][1] = -0.5; outputPntsMappedCoords[44][2] = -1.0;
      outputPntsMappedCoords[45][0] =  0.0; outputPntsMappedCoords[45][1] = -0.5; outputPntsMappedCoords[45][2] = -1.0;
      outputPntsMappedCoords[46][0] =  0.5; outputPntsMappedCoords[46][1] = -0.5; outputPntsMappedCoords[46][2] = -1.0;
      outputPntsMappedCoords[47][0] =  0.5; outputPntsMappedCoords[47][1] =  0.0; outputPntsMappedCoords[47][2] = -1.0;
      outputPntsMappedCoords[48][0] =  0.5; outputPntsMappedCoords[48][1] =  0.5; outputPntsMappedCoords[48][2] = -1.0;
      outputPntsMappedCoords[49][0] =  0.0; outputPntsMappedCoords[49][1] =  0.5; outputPntsMappedCoords[49][2] = -1.0;
      outputPntsMappedCoords[50][0] = -0.5; outputPntsMappedCoords[50][1] =  0.5; outputPntsMappedCoords[50][2] = -1.0;
      outputPntsMappedCoords[51][0] = -0.5; outputPntsMappedCoords[51][1] =  0.0; outputPntsMappedCoords[51][2] = -1.0;
      outputPntsMappedCoords[52][0] =  0.0; outputPntsMappedCoords[52][1] =  0.0; outputPntsMappedCoords[52][2] = -1.0;

      // ... face nodes 1
      outputPntsMappedCoords[53][0] = -0.5; outputPntsMappedCoords[53][1] = -1.0; outputPntsMappedCoords[53][2] = -0.5;
      outputPntsMappedCoords[54][0] =  0.0; outputPntsMappedCoords[54][1] = -1.0; outputPntsMappedCoords[54][2] = -0.5;
      outputPntsMappedCoords[55][0] =  0.5; outputPntsMappedCoords[55][1] = -1.0; outputPntsMappedCoords[55][2] = -0.5;
      outputPntsMappedCoords[56][0] =  0.5; outputPntsMappedCoords[56][1] = -1.0; outputPntsMappedCoords[56][2] =  0.0;
      outputPntsMappedCoords[57][0] =  0.5; outputPntsMappedCoords[57][1] = -1.0; outputPntsMappedCoords[57][2] =  0.5;
      outputPntsMappedCoords[58][0] =  0.0; outputPntsMappedCoords[58][1] = -1.0; outputPntsMappedCoords[58][2] =  0.5;
      outputPntsMappedCoords[59][0] = -0.5; outputPntsMappedCoords[59][1] = -1.0; outputPntsMappedCoords[59][2] =  0.5;
      outputPntsMappedCoords[60][0] = -0.5; outputPntsMappedCoords[60][1] = -1.0; outputPntsMappedCoords[60][2] =  0.0;
      outputPntsMappedCoords[61][0] =  0.0; outputPntsMappedCoords[61][1] = -1.0; outputPntsMappedCoords[61][2] =  0.0;
      // ... face nodes 2
      outputPntsMappedCoords[62][0] = 1.0; outputPntsMappedCoords[62][1] = -0.5; outputPntsMappedCoords[62][2] = -0.5;
      outputPntsMappedCoords[63][0] = 1.0; outputPntsMappedCoords[63][1] =  0.0; outputPntsMappedCoords[63][2] = -0.5;
      outputPntsMappedCoords[64][0] = 1.0; outputPntsMappedCoords[64][1] =  0.5; outputPntsMappedCoords[64][2] = -0.5;
      outputPntsMappedCoords[65][0] = 1.0; outputPntsMappedCoords[65][1] =  0.5; outputPntsMappedCoords[65][2] =  0.0;
      outputPntsMappedCoords[66][0] = 1.0; outputPntsMappedCoords[66][1] =  0.5; outputPntsMappedCoords[66][2] =  0.5;
      outputPntsMappedCoords[67][0] = 1.0; outputPntsMappedCoords[67][1] =  0.0; outputPntsMappedCoords[67][2] =  0.5;
      outputPntsMappedCoords[68][0] = 1.0; outputPntsMappedCoords[68][1] = -0.5; outputPntsMappedCoords[68][2] =  0.5;
      outputPntsMappedCoords[69][0] = 1.0; outputPntsMappedCoords[69][1] = -0.5; outputPntsMappedCoords[69][2] =  0.0;
      outputPntsMappedCoords[70][0] = 1.0; outputPntsMappedCoords[70][1] =  0.0; outputPntsMappedCoords[70][2] =  0.0;
      // ... face nodes 3
      outputPntsMappedCoords[71][0] =  0.5; outputPntsMappedCoords[71][1] = 1.0; outputPntsMappedCoords[71][2] = -0.5;
      outputPntsMappedCoords[72][0] =  0.0; outputPntsMappedCoords[72][1] = 1.0; outputPntsMappedCoords[72][2] = -0.5;
      outputPntsMappedCoords[73][0] = -0.5; outputPntsMappedCoords[73][1] = 1.0; outputPntsMappedCoords[73][2] = -0.5;
      outputPntsMappedCoords[74][0] = -0.5; outputPntsMappedCoords[74][1] = 1.0; outputPntsMappedCoords[74][2] =  0.0;
      outputPntsMappedCoords[75][0] = -0.5; outputPntsMappedCoords[75][1] = 1.0; outputPntsMappedCoords[75][2] =  0.5;
      outputPntsMappedCoords[76][0] =  0.0; outputPntsMappedCoords[76][1] = 1.0; outputPntsMappedCoords[76][2] =  0.5;
      outputPntsMappedCoords[77][0] =  0.5; outputPntsMappedCoords[77][1] = 1.0; outputPntsMappedCoords[77][2] =  0.5;
      outputPntsMappedCoords[78][0] =  0.5; outputPntsMappedCoords[78][1] = 1.0; outputPntsMappedCoords[78][2] =  0.0;
      outputPntsMappedCoords[79][0] =  0.0; outputPntsMappedCoords[79][1] = 1.0; outputPntsMappedCoords[79][2] =  0.0;
      // ... face nodes 4
      outputPntsMappedCoords[80][0] = -1.0; outputPntsMappedCoords[80][1] =  0.5; outputPntsMappedCoords[80][2] = -0.5;
      outputPntsMappedCoords[81][0] = -1.0; outputPntsMappedCoords[81][1] =  0.0; outputPntsMappedCoords[81][2] = -0.5;
      outputPntsMappedCoords[82][0] = -1.0; outputPntsMappedCoords[82][1] = -0.5; outputPntsMappedCoords[82][2] = -0.5;
      outputPntsMappedCoords[83][0] = -1.0; outputPntsMappedCoords[83][1] = -0.5; outputPntsMappedCoords[83][2] =  0.0;
      outputPntsMappedCoords[84][0] = -1.0; outputPntsMappedCoords[84][1] = -0.5; outputPntsMappedCoords[84][2] =  0.5;
      outputPntsMappedCoords[85][0] = -1.0; outputPntsMappedCoords[85][1] =  0.0; outputPntsMappedCoords[85][2] =  0.5;
      outputPntsMappedCoords[86][0] = -1.0; outputPntsMappedCoords[86][1] =  0.5; outputPntsMappedCoords[86][2] =  0.5;
      outputPntsMappedCoords[87][0] = -1.0; outputPntsMappedCoords[87][1] =  0.5; outputPntsMappedCoords[87][2] =  0.0;
      outputPntsMappedCoords[88][0] = -1.0; outputPntsMappedCoords[88][1] =  0.0; outputPntsMappedCoords[88][2] =  0.0;
      //---
      
      outputPntsMappedCoords[89][0] = -0.5; outputPntsMappedCoords[89][1] = -0.5; outputPntsMappedCoords[89][2] = 1.0;
      outputPntsMappedCoords[90][0] =  0.0; outputPntsMappedCoords[90][1] = -0.5; outputPntsMappedCoords[90][2] = 1.0;
      outputPntsMappedCoords[91][0] =  0.5; outputPntsMappedCoords[91][1] = -0.5; outputPntsMappedCoords[91][2] = 1.0;
      outputPntsMappedCoords[92][0] =  0.5; outputPntsMappedCoords[92][1] =  0.0; outputPntsMappedCoords[92][2] = 1.0;
      outputPntsMappedCoords[93][0] =  0.5; outputPntsMappedCoords[93][1] =  0.5; outputPntsMappedCoords[93][2] = 1.0;
      outputPntsMappedCoords[94][0] =  0.0; outputPntsMappedCoords[94][1] =  0.5; outputPntsMappedCoords[94][2] = 1.0;
      outputPntsMappedCoords[95][0] = -0.5; outputPntsMappedCoords[95][1] =  0.5; outputPntsMappedCoords[95][2] = 1.0;
      outputPntsMappedCoords[96][0] = -0.5; outputPntsMappedCoords[96][1] =  0.0; outputPntsMappedCoords[96][2] = 1.0;
      outputPntsMappedCoords[97][0] =  0.0; outputPntsMappedCoords[97][1] =  0.0; outputPntsMappedCoords[97][2] = 1.0;
      //---      
      
      outputPntsMappedCoords[98][0] = -0.5; outputPntsMappedCoords[98][1] = -0.5; outputPntsMappedCoords[98][2] = -0.5;
      outputPntsMappedCoords[99][0] =  0.0; outputPntsMappedCoords[99][1] = -0.5; outputPntsMappedCoords[99][2] = -0.5;
      outputPntsMappedCoords[100][0] =  0.5; outputPntsMappedCoords[100][1] = -0.5; outputPntsMappedCoords[100][2] = -0.5;
      outputPntsMappedCoords[101][0] =  0.5; outputPntsMappedCoords[101][1] =  0.0; outputPntsMappedCoords[101][2] = -0.5;
      outputPntsMappedCoords[102][0] =  0.5; outputPntsMappedCoords[102][1] =  0.5; outputPntsMappedCoords[102][2] = -0.5;
      outputPntsMappedCoords[103][0] =  0.0; outputPntsMappedCoords[103][1] =  0.5; outputPntsMappedCoords[103][2] = -0.5;
      outputPntsMappedCoords[104][0] = -0.5; outputPntsMappedCoords[104][1] =  0.5; outputPntsMappedCoords[104][2] = -0.5;
      outputPntsMappedCoords[105][0] = -0.5; outputPntsMappedCoords[105][1] =  0.0; outputPntsMappedCoords[105][2] = -0.5;
      outputPntsMappedCoords[106][0] =  0.0; outputPntsMappedCoords[106][1] =  0.0; outputPntsMappedCoords[106][2] = -0.5;
      //---          
        
      outputPntsMappedCoords[107][0] = -0.5; outputPntsMappedCoords[107][1] = -0.5; outputPntsMappedCoords[107][2] = 0.0;
      outputPntsMappedCoords[108][0] =  0.0; outputPntsMappedCoords[108][1] = -0.5; outputPntsMappedCoords[108][2] = 0.0;
      outputPntsMappedCoords[109][0] =  0.5; outputPntsMappedCoords[109][1] = -0.5; outputPntsMappedCoords[109][2] = 0.0;
      outputPntsMappedCoords[110][0] =  0.5; outputPntsMappedCoords[110][1] =  0.0; outputPntsMappedCoords[110][2] = 0.0;
      outputPntsMappedCoords[111][0] =  0.5; outputPntsMappedCoords[111][1] =  0.5; outputPntsMappedCoords[111][2] = 0.0;
      outputPntsMappedCoords[112][0] =  0.0; outputPntsMappedCoords[112][1] =  0.5; outputPntsMappedCoords[112][2] = 0.0;
      outputPntsMappedCoords[113][0] = -0.5; outputPntsMappedCoords[113][1] =  0.5; outputPntsMappedCoords[113][2] = 0.0;
      outputPntsMappedCoords[114][0] = -0.5; outputPntsMappedCoords[114][1] =  0.0; outputPntsMappedCoords[114][2] = 0.0;
      outputPntsMappedCoords[115][0] =  0.0; outputPntsMappedCoords[115][1] =  0.0; outputPntsMappedCoords[115][2] = 0.0;
      //---          
      
      outputPntsMappedCoords[116][0] = -0.5; outputPntsMappedCoords[116][1] = -0.5; outputPntsMappedCoords[116][2] = 0.5;
      outputPntsMappedCoords[117][0] =  0.0; outputPntsMappedCoords[117][1] = -0.5; outputPntsMappedCoords[117][2] = 0.5;
      outputPntsMappedCoords[118][0] =  0.5; outputPntsMappedCoords[118][1] = -0.5; outputPntsMappedCoords[118][2] = 0.5;
      outputPntsMappedCoords[119][0] =  0.5; outputPntsMappedCoords[119][1] =  0.0; outputPntsMappedCoords[119][2] = 0.5;
      outputPntsMappedCoords[120][0] =  0.5; outputPntsMappedCoords[120][1] =  0.5; outputPntsMappedCoords[120][2] = 0.5;
      outputPntsMappedCoords[121][0] =  0.0; outputPntsMappedCoords[121][1] =  0.5; outputPntsMappedCoords[121][2] = 0.5;
      outputPntsMappedCoords[122][0] = -0.5; outputPntsMappedCoords[122][1] =  0.5; outputPntsMappedCoords[122][2] = 0.5;
      outputPntsMappedCoords[123][0] = -0.5; outputPntsMappedCoords[123][1] =  0.0; outputPntsMappedCoords[123][2] = 0.5;
      outputPntsMappedCoords[124][0] =  0.0; outputPntsMappedCoords[124][1] =  0.0; outputPntsMappedCoords[124][2] = 0.5;
      break;
    }
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"ParCGNSHighOrderWriter::getOutputPntsMappedCoords Cell shape not supported");
    }
  }

  return outputPntsMappedCoords;
}
//////////////////////////////////////////////////////////////////////////////

void ParCGNSHighOrderWriter::writeP0CellCenteredSolution(
    MPI_Comm comm, int rank, int fileIndex, int baseIndex, 
    int index_zone, int index_sol, CFuint totNbCells, 
    CFuint globalNbCells, CFuint nbEqs,
    const std::vector<std::string>& varNames,
    const std::vector<std::string>& extraVarNames,
    const std::vector<std::string>& dh_varnames,
    Common::SafePtr<Framework::ConvectiveVarSet> updateVarSet,
    Common::SafePtr<Framework::DataHandleOutput> datahandle_output,
    Common::SafePtr<Framework::TopologicalRegionSet> trs,
    Framework::StdTrsGeoBuilder::GeoData& geoData,
    Common::SafePtr<std::vector<Framework::ElementTypeData>> elemType,
    CFuint nbrElemTypes)
{
  CFAUTOTRACE;
  
  // For P0, each cell has exactly one state at its center
  // No shape function evaluation needed - directly use cell state
  
  // Compute MPI ranges for cell-centered data
  CFuint localCellCount = totNbCells;
  CFuint startCell;
  
  MPI_Exscan(&localCellCount, &startCell, 1, MPI_UNSIGNED, MPI_SUM, comm);
  if (rank == 0) {
    startCell = 1;  // CGNS uses 1-based indexing
  } else {
    startCell += 1;
  }
  
  cgsize_t start_c = startCell;
  cgsize_t end_c = startCell + localCellCount - 1;
  
  // Allocate arrays for cell-centered data
  std::vector<std::vector<double>> cellDataArrays(
    nbEqs + extraVarNames.size() + dh_varnames.size(), 
    std::vector<double>(localCellCount)
  );
  
  // Helper states
  RealVector dimState(nbEqs);
  RealVector extraValues; // size will be set by VarSet if needed
  if (getMethodData().shouldPrintExtraValues() && !extraVarNames.empty()) {
    extraValues.resize(extraVarNames.size());
  }
  
  // Loop over element types
  CFuint cellIdx = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType) {
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    
    // Loop over cells in this element type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem) {
      geoData.idx = startIdx + iElem;
      GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();
      
      // Get cell states - for P0, there's exactly 1 state per cell
      vector<State*>* cellStates = cell->getStates();
      
      // Skip ghost cells
      bool isOwned = (!cellStates->empty() && (*cellStates)[0]->isParUpdatable());
      if (!isOwned) {
        m_stdTrsGeoBuilder.releaseGE();
        continue;
      }
      
      cf_assert(cellStates->size() == 1); // P0 has 1 state per cell
      const RealVector& cellState = *(*cellStates)[0];
      
      // Dimensionalize the state
      if (getMethodData().shouldPrintExtraValues()) {
        updateVarSet->setDimensionalValuesPlusExtraValues(cellState, dimState, extraValues);
      } else {
        updateVarSet->setDimensionalValues(cellState, dimState);
      }
      
      // Store solution variables
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        cellDataArrays[iEq][cellIdx] = dimState[iEq];
      }
      
      // Store extra variables
      for (CFuint iExtra = 0; iExtra < extraVarNames.size(); ++iExtra) {
        cellDataArrays[nbEqs + iExtra][cellIdx] = extraValues[iExtra];
      }
      
      // Store datahandle variables
      for (CFuint iVar = 0; iVar < dh_varnames.size(); ++iVar) {
        DataHandleOutput::DataHandleInfo var_info = datahandle_output->getStateData(iVar);
        CFuint var_var = var_info.first;
        CFuint var_nbvars = var_info.second;
        DataHandle<CFreal> var = var_info.third;
        
        cellDataArrays[nbEqs + extraVarNames.size() + iVar][cellIdx] = 
          var((*cellStates)[0]->getLocalID(), var_var, var_nbvars);
      }
      
      cellIdx++;
      m_stdTrsGeoBuilder.releaseGE();
    }
  }
  
  // Write cell-centered field data in parallel
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    if (cgp_field_write(fileIndex, baseIndex, index_zone, index_sol, 
                        RealDouble, varNames[iEq].c_str(), &fieldIndex) != CG_OK) {
      std::cerr << "Rank " << rank << " Error writing P0 field node for " 
                << varNames[iEq] << ": " << cg_get_error() << std::endl;
    }
    if (cgp_field_write_data(fileIndex, baseIndex, index_zone, index_sol, 
                             fieldIndex, &start_c, &end_c, cellDataArrays[iEq].data()) != CG_OK) {
      std::cerr << "Rank " << rank << " Error writing P0 field data for " 
                << varNames[iEq] << ": " << cg_get_error() << std::endl;
    }
  }
  
  // Write extra variables
  for (CFuint iExtra = 0; iExtra < extraVarNames.size(); ++iExtra) {
    if (cgp_field_write(fileIndex, baseIndex, index_zone, index_sol, 
                        RealDouble, extraVarNames[iExtra].c_str(), &fieldIndex) != CG_OK) {
      std::cerr << "Rank " << rank << " Error writing P0 extra field node for " 
                << extraVarNames[iExtra] << ": " << cg_get_error() << std::endl;
    }
    if (cgp_field_write_data(fileIndex, baseIndex, index_zone, index_sol, 
                             fieldIndex, &start_c, &end_c, 
                             cellDataArrays[nbEqs + iExtra].data()) != CG_OK) {
      std::cerr << "Rank " << rank << " Error writing P0 extra field data for " 
                << extraVarNames[iExtra] << ": " << cg_get_error() << std::endl;
    }
  }
  
  // Write datahandle variables
  for (CFuint iDh = 0; iDh < dh_varnames.size(); ++iDh) {
    if (cgp_field_write(fileIndex, baseIndex, index_zone, index_sol, 
                        RealDouble, dh_varnames[iDh].c_str(), &fieldIndex) != CG_OK) {
      std::cerr << "Rank " << rank << " Error writing P0 datahandle field node for " 
                << dh_varnames[iDh] << ": " << cg_get_error() << std::endl;
    }
    if (cgp_field_write_data(fileIndex, baseIndex, index_zone, index_sol, 
                             fieldIndex, &start_c, &end_c, 
                             cellDataArrays[nbEqs + extraVarNames.size() + iDh].data()) != CG_OK) {
      std::cerr << "Rank " << rank << " Error writing P0 datahandle field data for " 
                << dh_varnames[iDh] << ": " << cg_get_error() << std::endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////


