// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include <iomanip>

#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Common/OSystem.hh"

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"

#include "ParaViewWriter/ParaViewWriter.hh"
#include "ParaViewWriter/WriteSolutionHighOrder.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"
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

MethodCommandProvider<WriteSolutionHighOrder, ParaWriterData, ParaViewWriterModule>
writeSolutionHighOrderProvider("WriteSolutionHighOrder");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write ParaView file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionHighOrder::WriteSolutionHighOrder(const std::string& name) : ParaWriterCom(name)
{
  addConfigOptionsTo(this);

  m_fileFormatStr = "ASCII";
  setParameter("FileFormat",&m_fileFormatStr);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::execute()
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

const std::string WriteSolutionHighOrder::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::writeToBinaryFile()
{
 CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"Writing to binary files is not yet implemented in ParaViewWriter...");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  if (!getMethodData().onlySurface())
  {

  // get dimensionality, number of variables and reference length
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // variable that holds the indices of the variables that are vector components
  vector<CFuint> vectorComponentIdxs(0);
  /// @note this is a rather ugly piece of code, a check is made on the number of equations
  /// to avoid that vector components are searched when the physical model is a linear advection for instance
  /// this piece of code puts the velocity components (or the momentum components) in a vector
  /// the magnetic inductance vector B in the case of MHD is not put in a vector here like this!
  if (nbEqs >= 3)
  {
    switch (dim)
    {
      case DIM_2D:
      {
        vectorComponentIdxs.resize(2);
        vectorComponentIdxs[XX] = 1;
        vectorComponentIdxs[YY] = 2;
      } break;
      case DIM_3D:
      {
        vectorComponentIdxs.resize(3);
        vectorComponentIdxs[XX] = 1;
        vectorComponentIdxs[YY] = 2;
        vectorComponentIdxs[ZZ] = 3;
      } break;
      default:
      {
      }
    }
    cf_assert(dim == vectorComponentIdxs.size());
  }

  // variable that holds the indices of the scalar variables
  const CFuint nbVecComponents = vectorComponentIdxs.size();
  const CFuint nbScalars = nbEqs-nbVecComponents;
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

  // get ConvectiveVarSet
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  // get variable names
  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

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

  // loop over elements
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get element shape
    const CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();

    // get element solution order
    const CFuint solOrder = (*elemType)[iElemType].getSolOrder();

    // mapped coordinates of output points and cell-node connectivity
    const vector< RealVector >       outputPntsMappedCoords = getMethodData().getOutputPntsMappedCoords(shape,solOrder);
    const vector< vector< CFuint > > outputCellNodeConn     = getMethodData().getOutputCellNodeConn    (shape,solOrder);

    // number of output points
    const CFuint nbrOutPnts = outputPntsMappedCoords.size();

    // number of subcells
    const CFuint nbrSubCells = outputCellNodeConn.size();

    // get VTK cell type
    const CFuint vtkCellType = getMethodData().getVTKCellTypeID(shape,CFPolyOrder::ORDER1);

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
    GeometricEntity *const cell = geoBuilder->buildGE();
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
    geoBuilder->releaseGE();

    // coordinates and solutions in output nodes
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
      GeometricEntity *const cell = geoBuilder->buildGE();

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
          outputPntState[iPnt] += solShapeFuncs[iPnt][iState]*(*(*cellStates)[iState]);
        }
      }

      // open Piece element
      fout << "    <Piece NumberOfPoints=\"" << nbrOutPnts << "\" NumberOfCells=\"" << nbrSubCells << "\">\n";

      // open PointData element
//       fout << "      <PointData>\n";
      fout << "   <PointData Scalars=\"" << varNames[0] << "\">\n";

      // write the (velocity or momentum) vectors
      if (nbVecComponents > 0)
      {
        cf_assert(nbVecComponents >= 2);
        // open DataArray element
        fout << "        <DataArray Name=\"" << varNames[vectorComponentIdxs[1]] << "\" NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";
        fout << "          ";

        // loop over output points
        for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
        {
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

          // dimensionalize the state
          updateVarSet->setDimensionalValues(tempState, dimState);

          // write the current variable to the file
          for (CFuint iVecComp = 0; iVecComp < nbVecComponents; ++iVecComp)
          {
            fout << fixed << setprecision(12) << dimState[vectorComponentIdxs[iVecComp]] << " ";
          }
          for (CFuint iVecComp = nbVecComponents; iVecComp < 3; ++iVecComp)
          {
            fout << fixed << setprecision(1) << 0.0 << " ";
          }

        }
        fout << "\n";

        // close DataArray element
        fout << "        </DataArray>\n";
      }

      // loop over the scalars
      for (CFuint iScalar = 0; iScalar < nbScalars; ++iScalar)
      {
        // index of this scalar
        const CFuint iVar = scalarVarIdxs[iScalar];

        // open DataArray element
        fout << "        <DataArray Name=\"" << varNames[iVar] << "\" type=\"Float32\" format=\"ascii\">\n";
        fout << "          ";

        // loop over output points
        for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
        {
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

          // dimensionalize the state
          updateVarSet->setDimensionalValues(tempState, dimState);

          // write the current variable to the file
          fout << fixed << setprecision(12) << dimState[iVar] << " ";
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

          // loop over output points
          for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
          {
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
//             tempState.setSpaceCoordinates(&outputPntCoords[iPnt]);

            // dimensionalize the state
            updateVarSet->setDimensionalValuesPlusExtraValues(tempState, dimState, extraValues);

            // write the current variable to the file
            fout << fixed << setprecision(12) << extraValues[iVar] << " ";
          }
          fout << "\n";

          // close DataArray element
          fout << "        </DataArray>\n";
        }
      }

      // close PointData element
      fout << "      </PointData>\n";

      // open Points element
      fout << "      <Points>\n";

      // open DataArray element
      fout << "        <DataArray NumberOfComponents=\"" << 3 << "\" type=\"Float32\" format=\"ascii\">\n";
      fout << "          ";

      // loop over nodes to write coordinates
      for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
      {
        for (CFuint iCoor = 0; iCoor < dim; ++iCoor)
        {
          fout << fixed << setprecision(12) << outputPntCoords[iPnt][iCoor]*refL << " ";
        }
        for (CFuint iCoor = dim; iCoor < 3; ++iCoor)
        {
          fout << fixed << setprecision(1) << 0.0 << " ";
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

      // loop over cells to write cell-node connectivity
      for (CFuint iCell = 0; iCell < nbrSubCells; ++iCell)
      {
        const CFuint nbrNodes = outputCellNodeConn[iCell].size();
        // node ordering for one cell is the same for VTK as in COOLFluiD
        for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
        {
          fout << outputCellNodeConn[iCell][iNode] << " ";
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
      for (CFuint iCell = 0; iCell < nbrSubCells; ++iCell)
      {
        cellEndOffSet += outputCellNodeConn[iCell].size();
        fout << cellEndOffSet << " ";
      }
      fout << "\n";

      // close DataArray element (for offsets in cell-node connectivity)
      fout << "        </DataArray>\n";

      // open DataArray element (for cell types)
      fout << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
      fout << "          ";

      // loop over element types to write cell types
      for (CFuint iCell = 0; iCell < nbrSubCells; ++iCell)
      {
        fout << vtkCellType << " ";
      }
      fout << "\n";

      // close DataArray element (for cell types)
      fout << "        </DataArray>\n";

      // close Cells element
      fout << "      </Cells>\n";

      // close Piece element
      fout << "    </Piece >\n";


      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }

  // close UnstructuredGrid element
  fout << "  </UnstructuredGrid>\n";

  // close VTKFile element
  fout << "</VTKFile>\n";

  // close the file
  fout.close();

  } // if only surface

  // write boundary surface data
  writeBoundarySurface();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHighOrder::writeBoundarySurface()
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

void WriteSolutionHighOrder::setup()
{
  CFAUTOTRACE;

  ParaWriterCom::setup();

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  updateVarSet->setup();

  SafePtr< GeometricEntityPool< StdTrsGeoBuilder > > cellBuilder = getMethodData().getStdTrsGeoBuilder();

  cellBuilder->setup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionHighOrder::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD
