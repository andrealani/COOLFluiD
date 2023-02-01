// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include <iomanip>

#include "Common/COOLFluiD.hh"
#ifdef CF_HAVE_UNISTD_H
  extern "C"
  {
    #include "TecplotWriter/TECXXX.h"
  }
#endif

#include "Common/PE.hh"

#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Common/OSystem.hh"

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
// #include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/StdTrsGeoBuilder.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/WriteSolutionHighOrder.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataHandleOutput.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteSolutionHighOrder, TecWriterData, TecplotWriterModule>
writeSolutionHighOrderProvider("WriteSolutionHighOrder");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionHighOrder::WriteSolutionHighOrder(const std::string& name) : TecWriterCom(name),
    m_stdTrsGeoBuilder()
{
  addConfigOptionsTo(this);

  _fileFormatStr = "ASCII";
  setParameter("FileFormat",&_fileFormatStr);

}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::execute()
{
  CFLog(INFO, "Writing solution to: " << getMethodData().getFilename().string() << "\n");

  if(_fileFormatStr == "ASCII"){
    writeToFile(getMethodData().getFilename());
  }
  else
  {
    cf_assert(_fileFormatStr == "BINARY");

    ///@todo change this to use the tecplot library
    ///this is slow and NOT portable but at least, it takes less space
    writeToFile("tmp");
    std::string transformFile = "$TECHOME/bin/preplot tmp " + getMethodData().getFilename().string();
    CFLog(INFO, transformFile << "\n");
    
    OSystem::getInstance().executeCommand(transformFile);

//     writeToBinaryFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

const std::string WriteSolutionHighOrder::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::writeToBinaryFile()
{
 CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

//   SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if (!getMethodData().onlySurface())
  {
  // get dimensionality, number of equations and reference length
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // write Tecplot header
  fout << "TITLE      =  \"Unstructured grid data\"" << "\n";
  fout << "VARIABLES  = ";
  for (CFuint i = 0; i < dim; ++i) {
    fout << " \"x" << i << '\"';
  }

  // get convective variable set
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  // get variable names
  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  // write variable names
  for (CFuint i = 0 ;  i < nbEqs; ++i) {
    std::string n = varNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    fout << " " << n;
  }
  if (getMethodData().shouldPrintExtraValues()) {
    vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
    for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
      fout << " " << extraVarNames[i];
    }
  }
  
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();
  std::vector< std::string > dh_varnames = datahandle_output->getVarNames();
  for (CFuint iVar = 0; iVar < dh_varnames.size(); ++iVar)
  {
    std::string n = dh_varnames[iVar];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    fout << " " << n;
  }

  fout << "\n";

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");

  // prepares to loop over cells by getting the GeometricEntityPool
  StdTrsGeoBuilder::GeoData& geoData = m_stdTrsGeoBuilder.getDataGE();
  geoData.trs = trs;

  // loop over elements
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get element shape
    const CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();

    // get element solution order
    const CFuint solOrder = (*elemType)[iElemType].getSolOrder();

    // mapped coordinates of output points and cell-node connectivity
    const vector< RealVector >       outputPntsMappedCoords = getOutputPntsMappedCoords(shape,solOrder);
    const vector< vector< CFuint > > outputCellNodeConn     = getOutputCellNodeConn    (shape,solOrder);

    // number of output points
    const CFuint nbrOutPnts = outputPntsMappedCoords.size();

    // number of subcells
    const CFuint nbrSubCells = outputCellNodeConn.size();

    // number of states in element type
    const CFuint nbrStates = (*elemType)[iElemType].getNbStates();

    // number of nodes in element type
    const CFuint nbrNodes = (*elemType)[iElemType].getNbNodes();

    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // evaluate the basis functions in the output points
    geoData.idx = cellIdx;
    GeometricEntity *const cell = m_stdTrsGeoBuilder.buildGE();
    vector< RealVector > solShapeFuncs;
    vector< RealVector > geoShapeFuncs;
    for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
    {
      RealVector solShapeFunc =
          cell->computeShapeFunctionAtMappedCoord   (outputPntsMappedCoords[iPnt]);
      solShapeFuncs.push_back(solShapeFunc);
      RealVector geoShapeFunc =
          cell->computeGeoShapeFunctionAtMappedCoord(outputPntsMappedCoords[iPnt]);
      geoShapeFuncs.push_back(geoShapeFunc);
    }
    //release the GeometricEntity
    m_stdTrsGeoBuilder.releaseGE();

    // variable for coordinates and solutions in output nodes
    vector< Node  > outputPntCoords;
    vector< State > outputPntState;
    for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
    {
      RealVector aux1(dim);
      Node  coords(aux1,false);
      outputPntCoords.push_back(coords);
      RealVector aux2(nbEqs);
      State state(aux2,false);
      outputPntState.push_back(state);
    }

    // some helper states
    RealVector dimState(nbEqs);
    RealVector extraValues; // size will be set in the VarSet
    State tempState;

      // loop over elements
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

      // evaluate node coordinates and states at the output points
      for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
      {
        outputPntCoords[iPnt] = 0.0;
        for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
        {
          outputPntCoords[iPnt] += geoShapeFuncs[iPnt][iNode]*(*(*cellNodes)[iNode]);
        }

        outputPntState[iPnt]  = 0.0;
        for (CFuint iState = 0; iState < nbrStates; ++iState)
	{  
	  // CFLog(INFO, "state[" << iState << "] = " << (*(*cellStates)[iState]) << ", solShapeFuncs[" << iPnt << "] = " << solShapeFuncs[iPnt][iState] <<"\n");
	  outputPntState[iPnt] += solShapeFuncs[iPnt][iState]*(*(*cellStates)[iState]);
        }
      }

      // create new zone for this cell
      fout << "ZONE "
          << "  T=\"ZONE " << cellIdx << "\""
          << ", N=" << nbrOutPnts
          << ", E=" << nbrSubCells
          << ", F=FEPOINT"
          << ", ET=" << getTecplotCellShape(shape,1);
//          if (getMethodData().getAppendAuxData())
//            fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank() << "\""
//                 << ", AUXDATA TRS=\"" << trs->getName() << "\""
//                 << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
//                 << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
//                 << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
//                 << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
//                 << flush;
      fout <<  "\n" << flush;

      // write node coordinates and states
      for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
      {
        // write coordinates
        for (CFuint iDim = 0; iDim < dim; ++iDim) {
          fout << setw(20) << scientific << setprecision(12)
               << outputPntCoords[iPnt][iDim]*refL << " ";
        }

        // get state in this point
        const RealVector& nodalState = outputPntState[iPnt];

        // copy to temporary state variable
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
	{
	  tempState[iEq] = nodalState[iEq];
        }
		
        // set temporary state ID
        CFLogDebugMin("temporary state ID not set in WriteSolutionHighOrder!!");
//           const CFuint stateID = nodalStates.getStateLocalID(iNode);
//           tempState.setLocalID(stateID);

        // set the node in the temporary state
        CFLogDebugMin("temporary state node not set in WriteSolutionHighOrder!!");
//           tempState.setSpaceCoordinates(&outputPntCoords[iPnt]);


        // dimensionalize the solution and compute extra values
        if (getMethodData().shouldPrintExtraValues())
        {
          updateVarSet->setDimensionalValuesPlusExtraValues
              (tempState, dimState, extraValues);
        }
        else
        {
          updateVarSet->setDimensionalValues(tempState, dimState);
        }

        // write state
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          fout << setw(20) << scientific << setprecision(12)
               << dimState[iEq] << " ";
        }

        // write extra values
        const CFuint nbrExtraVals = extraValues.size();
        for (CFuint iVal = 0; iVal < nbrExtraVals; ++iVal)
        {
          fout << setw(20) << scientific << setprecision(12)
               << extraValues[iVal] << " ";
        }
	
	      // print datahandles with state based data
  {

    // loop over the state based variables

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
	
      // write the current variable to the file
      fout << setw(20) << scientific << setprecision(12)
             << outputPntStateSockets[iPnt] << " ";      
    }
  }
  fout << "\n";
      }
      fout << "\n";

      // write connectivity
      for (CFuint iCell = 0; iCell < nbrSubCells; ++iCell)
      {
        const CFuint nbrNodes = outputCellNodeConn[iCell].size();
        // node ordering for one cell is the same for VTK as in COOLFluiD
        for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
        {
          fout << outputCellNodeConn[iCell][iNode]+1 << " ";
        }
        fout << "\n";
      }
      fout << "\n";

      //release the GeometricEntity
      m_stdTrsGeoBuilder.releaseGE();
    }
  }

  fout.close();

  } // if only surface

  // write boundary surface data
  writeBoundarySurface();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();

  // setup geobuilder
  m_stdTrsGeoBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::writeBoundarySurface()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionHighOrder::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::string WriteSolutionHighOrder::getTecplotCellShape(CFGeoShape::Type shape,CFuint geoOrder)
{
  switch (shape)
  {
    case CFGeoShape::TRIAG:
    {
      return "TRIANGLE";
    }
    case CFGeoShape::QUAD:
    {
      return "QUADRILATERAL";
    }
    case CFGeoShape::TETRA:
    {
      return "TETRAHEDRON";
    }
    case CFGeoShape::PYRAM: // same as hexahedra but with repeated nodes
    case CFGeoShape::PRISM: // same as hexahedra but with repeated nodes
    case CFGeoShape::HEXA:
    {
      return "BRICK";
    }
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Cell shape not supported by tecplot (20-05-2008)");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > WriteSolutionHighOrder::getOutputPntsMappedCoords(CFGeoShape::Type shape,CFuint solOrder)
{
  // output variable
  vector< RealVector > nodeMappedCoords;

  // for zeroth order solution polynomial, use same points as for first order
  if (solOrder == 0)
  {
    solOrder = 1;
  }

  // number of points needed for representing a polynomial of order degree solOrder
  const CFuint nbrNodes1D = 1*(solOrder + 1);

  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        RealVector coords(1);
        coords[KSI] = -1.0 + iKsi*2.0/solOrder;
        nodeMappedCoords.push_back(coords);
      }
    } break;
    case CFGeoShape::TRIAG:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal ksi = iKsi*1.0/solOrder;
        for (CFuint iEta = 0; iEta < nbrNodes1D-iKsi; ++iEta)
        {
          RealVector coords(2);
          coords[KSI] = ksi;
          coords[ETA] = iEta*1.0/solOrder;
          nodeMappedCoords.push_back(coords);
        }
      }
    } break;
    case CFGeoShape::QUAD:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal ksi = -1.0 + iKsi*2.0/(nbrNodes1D-1);
        for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
        {
          RealVector coords(2);
          coords[KSI] = ksi;
          coords[ETA] = -1.0 + iEta*2.0/(nbrNodes1D-1);
          nodeMappedCoords.push_back(coords);
        }
      }
    } break;
    case CFGeoShape::TETRA:
    {
      /// @warn: not very efficient
      std::vector<std::vector<CFuint> > cellsNodesConn;

      RealVector coords(3);
      coords[KSI] = 0.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 1.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 0.0;
      coords[ETA] = 1.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 0.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 1.0;
      nodeMappedCoords.push_back(coords);

      std::vector<CFuint> nodes_index(nodeMappedCoords.size());
      for (CFuint i = 0; i < nodes_index.size(); ++i) 
      {
        nodes_index[i] = i;
      }
      vector<vector<CFuint> > cells;
      vector<CFuint> cell(4);
      cell[0] = 0;
      cell[1] = 1;
      cell[2] = 2;
      cell[3] = 3;
      cells.push_back(cell);

      for (int i = 0; i < solOrder-1; i++) 
      {
        std::vector<std::vector<CFuint> > new_cells;
        for (int j = 0; j < cells.size(); j++) 
        {
          vector<RealVector> nodes(4);
          for (CFuint c = 0; c < 4; c++)
          {
            nodes[c].resize(3); 
          }
          nodes[0] = nodeMappedCoords[cells[j][0]];
          nodes[1] = nodeMappedCoords[cells[j][1]];
          nodes[2] = nodeMappedCoords[cells[j][2]];
          nodes[3] = nodeMappedCoords[cells[j][3]];

          vector<RealVector> midpoints = compute_midpoints(nodes);

          vector<vector<RealVector> > sub_tetrahedra = subdivide_tetrahedron(nodes, midpoints);

          for (int k = 0; k < sub_tetrahedra.size(); k++) 
          {
            for (int l = 0; l < sub_tetrahedra[k].size(); l++) 
            {
              CFuint node_index = std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][l]) - nodeMappedCoords.begin();
              if (node_index == nodeMappedCoords.size()) 
              {
                  nodeMappedCoords.push_back(sub_tetrahedra[k][l]);
                  node_index = nodeMappedCoords.size() - 1;
              }
            }
            vector<CFuint> new_cell(4);
            new_cell[0] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][0]) - nodeMappedCoords.begin());
            new_cell[1] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][1]) - nodeMappedCoords.begin());
            new_cell[2] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][2]) - nodeMappedCoords.begin());
            new_cell[3] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][3]) - nodeMappedCoords.begin());
            new_cells.push_back(new_cell);
          }
        }
        cells = new_cells;
      }

      cellsNodesConn = cells;


    } break;
    case CFGeoShape::PYRAM:
    {
      throw Common::NotImplementedException (FromHere(),"WriteSolutionHighOrder::getOutputPntsMappedCoords() for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"WriteSolutionHighOrder::getOutputPntsMappedCoords() for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal ksi = -1.0 + iKsi*2.0/solOrder;
        for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
        {
          const CFreal eta = -1.0 + iEta*2.0/solOrder;
          for (CFuint iZta = 0; iZta < nbrNodes1D; ++iZta)
          {
            RealVector coords(3);
            coords[KSI] = ksi;
            coords[ETA] = eta;
            coords[ZTA] = -1.0 + iZta*2.0/solOrder;
            nodeMappedCoords.push_back(coords);
          }
        }
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
    }
  }

  return nodeMappedCoords;
}

//////////////////////////////////////////////////////////////////////////////

vector< vector< CFuint > > WriteSolutionHighOrder::getOutputCellNodeConn(CFGeoShape::Type shape,CFuint solOrder)
{
  // output variable
  vector< vector< CFuint > > cellsNodesConn;

  // for zeroth order solution polynomial, use same points as for first order
  if (solOrder == 0)
  {
    solOrder = 1;
  }

  // number of points needed for representing a polynomial of order degree solOrder
  const CFuint nbrNodes1D = 1*(solOrder + 1);

  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        vector< CFuint > cellNodesConn(2);
        cellNodesConn[0] = iKsi;
        cellNodesConn[1] = iKsi+1;
        cellsNodesConn.push_back(cellNodesConn);
      }
    } break;
    case CFGeoShape::TRIAG:
    {
      CFuint nodeIdx = 0;
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < solOrder-iKsi-1; ++iEta)
        {
          vector< CFuint > cellNodesConn(3);
          cellNodesConn[0] = nodeIdx + iEta;
          cellNodesConn[1] = nodeIdx + nbrNodes1D - iKsi + iEta;
          cellNodesConn[2] = nodeIdx + iEta + 1;
          cellsNodesConn.push_back(cellNodesConn);

          cellNodesConn[0] = nodeIdx + nbrNodes1D - iKsi + iEta;
          cellNodesConn[1] = nodeIdx + nbrNodes1D - iKsi + iEta + 1;
          cellNodesConn[2] = nodeIdx + iEta + 1;
          cellsNodesConn.push_back(cellNodesConn);
        }
        
        vector< CFuint > cellNodesConn(3);
        cellNodesConn[0] = nodeIdx + solOrder-iKsi-1;
        cellNodesConn[1] = nodeIdx + nbrNodes1D -iKsi + solOrder - 1 - iKsi;
        cellNodesConn[2] = nodeIdx + solOrder-iKsi;
        cellsNodesConn.push_back(cellNodesConn);
        
        nodeIdx += nbrNodes1D-iKsi;
      }
    } break;
    case CFGeoShape::QUAD:
    {
      for (CFuint iKsi = 0; iKsi < (nbrNodes1D-1); ++iKsi)
      {
        for (CFuint iEta = 0; iEta < (nbrNodes1D-1); ++iEta)
        {
          vector< CFuint > cellNodesConn(4);
          cellNodesConn[0] = (iKsi  )*nbrNodes1D + iEta  ;
          cellNodesConn[1] = (iKsi+1)*nbrNodes1D + iEta  ;
          cellNodesConn[2] = (iKsi+1)*nbrNodes1D + iEta+1;
          cellNodesConn[3] = (iKsi  )*nbrNodes1D + iEta+1;
          cellsNodesConn.push_back(cellNodesConn);
        }
      }
    } break;
    case CFGeoShape::TETRA:
    {
      /// @warn: not very effcient

      std::vector<RealVector> nodeMappedCoords;
      RealVector coords(3);
      coords[KSI] = 0.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 1.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 0.0;
      coords[ETA] = 1.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 0.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 1.0;
      nodeMappedCoords.push_back(coords);
      std::vector<CFuint> nodes_index(nodeMappedCoords.size());
      for (CFuint i = 0; i < nodes_index.size(); ++i) 
      {
        nodes_index[i] = i;
      }
      vector<vector<CFuint> > cells;
      vector<CFuint> cell(4);
      cell[0] = 0;
      cell[1] = 1;
      cell[2] = 2;
      cell[3] = 3;
      cells.push_back(cell);

      for (int i = 0; i < solOrder-1; i++) 
      {
        std::vector<std::vector<CFuint> > new_cells;
        for (int j = 0; j < cells.size(); j++) 
        {
          vector<RealVector> nodes(4);
          for (CFuint c = 0; c < 4; c++)
          {
            nodes[c].resize(3); 
          }
          nodes[0] = nodeMappedCoords[cells[j][0]];
          nodes[1] = nodeMappedCoords[cells[j][1]];
          nodes[2] = nodeMappedCoords[cells[j][2]];
          nodes[3] = nodeMappedCoords[cells[j][3]];
          vector<RealVector> midpoints = compute_midpoints(nodes);

          vector<vector<RealVector> > sub_tetrahedra = subdivide_tetrahedron(nodes, midpoints);

          for (int k = 0; k < sub_tetrahedra.size(); k++) 
          {
            for (int l = 0; l < sub_tetrahedra[k].size(); l++) 
            {
              CFuint node_index = std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][l]) - nodeMappedCoords.begin();
              if (node_index == nodeMappedCoords.size()) 
              {
                  nodeMappedCoords.push_back(sub_tetrahedra[k][l]);
                  node_index = nodeMappedCoords.size() - 1;
              }
            }
            vector<CFuint> new_cell(4);
            new_cell[0] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][0]) - nodeMappedCoords.begin());
            new_cell[1] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][1]) - nodeMappedCoords.begin());
            new_cell[2] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][2]) - nodeMappedCoords.begin());
            new_cell[3] = (std::find(nodeMappedCoords.begin(), nodeMappedCoords.end(), sub_tetrahedra[k][3]) - nodeMappedCoords.begin());
            new_cells.push_back(new_cell);          
          }
        }
          cells = new_cells;
      }

      cellsNodesConn = cells;

    } break;
    case CFGeoShape::PYRAM:
    {
      throw Common::NotImplementedException (FromHere(),"WriteSolutionHighOrder::getOutputCellNodeConn() for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"WriteSolutionHighOrder::getOutputCellNodeConn() for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      const CFuint nbrNodes1DSq = nbrNodes1D*nbrNodes1D;
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < solOrder; ++iEta)
        {
          for (CFuint iZta = 0; iZta < solOrder; ++iZta)
          {
            vector< CFuint > cellNodesConn(8);
            cellNodesConn[0] = (iKsi  )*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta  ;
            cellNodesConn[1] = (iKsi+1)*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta  ;
            cellNodesConn[2] = (iKsi+1)*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta  ;
            cellNodesConn[3] = (iKsi  )*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta  ;
            cellNodesConn[4] = (iKsi  )*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta+1;
            cellNodesConn[5] = (iKsi+1)*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta+1;
            cellNodesConn[6] = (iKsi+1)*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta+1;
            cellNodesConn[7] = (iKsi  )*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta+1;
            cellsNodesConn.push_back(cellNodesConn);
          }
        }
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
    }
  }

  return cellsNodesConn;
}

//////////////////////////////////////////////////////////////////////////////

vector<RealVector> compute_midpoints(vector<RealVector> tetrahedron) 
{
  vector<RealVector> midpoints;
  for (CFuint i = 0; i < 4; i++) 
  {
    for (CFuint j = i + 1; j < 4; j++) 
    {
      RealVector midpoint(3);
      for (CFuint k = 0; k < 3; k++) 
      {
        midpoint[k] = (tetrahedron[i][k] + tetrahedron[j][k]) / 2;
      }
      midpoints.push_back(midpoint);
    }
  }
  return midpoints;
}

//////////////////////////////////////////////////////////////////////////////

vector<vector<RealVector> > subdivide_tetrahedron(vector<RealVector> tetrahedron, vector<RealVector> midpoints) 
{
  vector<vector<RealVector> > sub_tetrahedra;
  vector<RealVector> sub_tetrahedron1(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron1[c].resize(3); 
  }
  sub_tetrahedron1[0] = tetrahedron[0];
  sub_tetrahedron1[1] = midpoints[0];
  sub_tetrahedron1[2] = midpoints[1];
  sub_tetrahedron1[3] = midpoints[2];
  vector<RealVector> sub_tetrahedron2(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron2[c].resize(3); 
  }
  sub_tetrahedron2[0] = tetrahedron[1];
  sub_tetrahedron2[1] = midpoints[3];
  sub_tetrahedron2[2] = midpoints[0];
  sub_tetrahedron2[3] = midpoints[4];
  vector<RealVector> sub_tetrahedron3(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron3[c].resize(3); 
  }
  sub_tetrahedron3[0] = tetrahedron[2];
  sub_tetrahedron3[1] = midpoints[1];
  sub_tetrahedron3[2] = midpoints[3];
  sub_tetrahedron3[3] = midpoints[5];
  vector<RealVector> sub_tetrahedron4(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron4[c].resize(3); 
  }
  sub_tetrahedron4[0] = tetrahedron[3];
  sub_tetrahedron4[1] = midpoints[4];
  sub_tetrahedron4[2] = midpoints[5];
  sub_tetrahedron4[3] = midpoints[2];
  vector<RealVector> sub_tetrahedron5(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron5[c].resize(3); 
  }
  sub_tetrahedron5[0] = midpoints[0];
  sub_tetrahedron5[1] = midpoints[1];
  sub_tetrahedron5[2] = midpoints[3];
  sub_tetrahedron5[3] = midpoints[4];
  vector<RealVector> sub_tetrahedron6(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron6[c].resize(3); 
  }
  sub_tetrahedron6[0] = midpoints[0];
  sub_tetrahedron6[1] = midpoints[1];
  sub_tetrahedron6[2] = midpoints[4];
  sub_tetrahedron6[3] = midpoints[2];
  vector<RealVector> sub_tetrahedron7(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron7[c].resize(3); 
  }
  sub_tetrahedron7[0] = midpoints[2];
  sub_tetrahedron7[1] = midpoints[5];
  sub_tetrahedron7[2] = midpoints[4];
  sub_tetrahedron7[3] = midpoints[1];
  vector<RealVector> sub_tetrahedron8(4);
  for (CFuint c = 0; c < 4; c++)
  {
    sub_tetrahedron8[c].resize(3); 
  }
  sub_tetrahedron8[0] = midpoints[1];
  sub_tetrahedron8[1] = midpoints[5];
  sub_tetrahedron8[2] = midpoints[4];
  sub_tetrahedron8[3] = midpoints[3];
  sub_tetrahedra.push_back(sub_tetrahedron1);
  sub_tetrahedra.push_back(sub_tetrahedron2);
  sub_tetrahedra.push_back(sub_tetrahedron3);
  sub_tetrahedra.push_back(sub_tetrahedron4);
  sub_tetrahedra.push_back(sub_tetrahedron5);
  sub_tetrahedra.push_back(sub_tetrahedron6);
  sub_tetrahedra.push_back(sub_tetrahedron7);
  sub_tetrahedra.push_back(sub_tetrahedron8);
  return sub_tetrahedra;
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
