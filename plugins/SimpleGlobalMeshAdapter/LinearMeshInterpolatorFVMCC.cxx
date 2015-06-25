// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "LinearMeshInterpolatorFVMCC.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/PhysicalModel.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LinearMeshInterpolatorFVMCC, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> linearMeshInterpolatorFVMCCProvider("LinearMeshInterpolatorFVMCC");

//////////////////////////////////////////////////////////////////////////////

LinearMeshInterpolatorFVMCC::LinearMeshInterpolatorFVMCC(const std::string& name)  :
  LinearMeshInterpolator(name),
  socket_otherNstates("nstates"),
  socket_otherGstates("gstates"),
  socket_otherNodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolatorFVMCC::configure ( Config::ConfigArgs& args )
{
  LinearMeshInterpolator::configure(args);

  socket_otherNstates.setDataSocketNamespace(getMethodData().getOtherNamespace());
  socket_otherGstates.setDataSocketNamespace(getMethodData().getOtherNamespace());
  socket_otherNodes.setDataSocketNamespace(getMethodData().getOtherNamespace());
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolatorFVMCC::setup()
{
  CFAUTOTRACE;

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  cf_assert(spaceMethod.isNotNull());

  LinearMeshInterpolator::setup();

}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolatorFVMCC::execute()
{
  CFAUTOTRACE;

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  cf_assert(spaceMethod.isNotNull());

  spaceMethod->extrapolateStatesToNodes();

  LinearMeshInterpolator::execute();

}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolatorFVMCC::BoxClassification()
{
  CFAUTOTRACE;

  //Create a structured mesh dividing the domain into rectangular zones
  //which contains only a subset of cells
  //You do this once and it is then faster to search

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  std::string otherNamespace = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(otherNamespace);
  
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  GeometricEntityPool<CellTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setupInNamespace(otherNamespace);

  SafePtr<CellTrsGeoBuilder> geoBuilderCellPtr = geoBuilderCell.getGeoBuilder();
  geoBuilderCellPtr->setDataSockets(socket_otherStates, socket_otherGstates, socket_otherNodes);

  CellTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();
  geoDataCell.trs = cells;

  const CFuint maxNbNodesInCell =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->Statistics().getMaxNbNodesInCell();
  std::vector<CFuint> cellBoxes(maxNbNodesInCell);

  //Loop over the cells in the whole domain and put the cells in the boxes
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    geoDataCell.idx = iCell;
    GeometricEntity& cell = *geoBuilderCell.buildGE();
    const CFuint nbNodesInCell = cell.nbNodes();

//   CFout << "Cell: " << iCell <<"\n";

    //initialize to extreme values
    for(CFuint iDim = 0; iDim < nbDim; iDim++)
    {
      _minCoord[iDim] = MathTools::MathConsts::CFrealMax();
      _maxCoord[iDim] = -MathTools::MathConsts::CFrealMax();
    }

    //Compute the max and min coordinates
    for(CFuint iNode=0; iNode < nbNodesInCell; iNode++)
    {
      _coord = *(cell.getNode(iNode));
//  CFout << "Coord: " << _coord <<"\n";
      for(CFuint iDim = 0; iDim < nbDim; iDim++)
      {
        _minCoord[iDim] = std::min(_minCoord[iDim], _coord[iDim]);
        _maxCoord[iDim] = std::max(_maxCoord[iDim], _coord[iDim]);
      }
    }

//  CFout << "Min Coord: " << _minCoord <<"\n";
//  CFout << "Max Coord: " << _maxCoord <<"\n";
//  CFout << "Trying to find the boxes that contain the cell\n";

    findBoxesFromMinMax();

    for(CFuint iBox=0;iBox<_listBoxes.size(); ++iBox)
    {
      _boxes[_listBoxes[iBox]].push_back(iCell);
//   CFout << "Put in Box: " << _listBoxes[iBox] << "\n";
    }

    //release the GeometricEntity
    geoBuilderCell.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolatorFVMCC::SearchInBoxes()
{

  CFAUTOTRACE;

  DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  DataHandle<RealVector > otherNstates = socket_otherNstates.getDataHandle();
  DataHandle<Framework::State*, Framework::GLOBAL> otherStates = socket_otherStates.getDataHandle();

// CFout << "Searching for state at coord: " << _coord <<"\n";
  const CFuint boxIndex = findBoxFromCoord();
  const CFuint nbCellsInBox = _boxes[boxIndex].size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->getTrs("InnerCells");

  GeometricEntityPool<CellTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setupInNamespace(name);

  SafePtr<CellTrsGeoBuilder> geoBuilderCellPtr = geoBuilderCell.getGeoBuilder();
  geoBuilderCellPtr->setDataSockets(socket_otherStates, socket_otherGstates, socket_otherNodes);

  CellTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();
  geoDataCell.trs = cells;

  //Loop over the cells in the box
  bool elementFound(false);
  for (CFuint iCell = 0; (iCell < nbCellsInBox) && (!elementFound); ++iCell) {

    // build the GeometricEntity
    const CFuint cellID = _boxes[boxIndex][iCell];
    geoDataCell.idx = cellID;
    GeometricEntity *const currCell = geoBuilderCell.buildGE();

    //check if coord is in cell
    elementFound = currCell->isInElement(_coord);
    if(elementFound)
    {
      _shapeFunctions = currCell->computeGeoShapeFunctionAtCoord(_coord);

      std::vector<Node*>* cellNodes = currCell->getNodes();
//      std::vector<State*>* cellStates = currCell->getStates();
      const CFuint nbNodes = currCell->nbNodes();
      const CFuint nbStates = currCell->nbStates();

      //this is for FVM only!!!
      cf_assert(nbStates == 1);

//        for (CFuint j = 0; j < nbEqs; ++j)
//        {
//          _interpolatedState[j] = (*((*cellStates)[0]))[j];
//        }


      for (CFuint j = 0; j < nbEqs; ++j)
      {
        CFuint nodeID = ((*cellNodes)[0])->getLocalID();
        _interpolatedState[j] = _shapeFunctions[0] * (otherNstates[nodeID])[j];
        for (CFuint k = 1; k < nbNodes; ++k)
        {
          nodeID = ((*cellNodes)[k])->getLocalID();
          _interpolatedState[j] += _shapeFunctions[k] * (otherNstates[nodeID])[j];
        }
      }

      (states[_stateID])->copyData(_interpolatedState);
    }

    //release the GeometricEntity
    geoBuilderCell.releaseGE();
  }

  //if false: compute the closest point (projection)
  if(!elementFound){
//     CFout << "No element was found for state: " <<_coord <<" in Box: "<< boxIndex<<"\n";
//     CFout << "It must be outside of the old domain...\n";
//     CFout << "Let's interpolate using the nearest nodes...\n";

    for (CFuint iCell = 0; (iCell < nbCellsInBox); ++iCell) {
//      CFout << "Containing cell: " << _boxes[boxIndex][iCell] <<"\n";
    }
    computeClosestValue();
    (states[_stateID])->copyData(_interpolatedState);
  }

}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolatorFVMCC::computeClosestValue()
{
//on which boundary should we look??
//Look on the boundaries with the same
//cf_assert(false);

  DataHandle<RealVector > otherStates = socket_otherNstates.getDataHandle();
  DataHandle<Node*,GLOBAL> otherNodes = socket_otherNodes.getDataHandle();
///@todo compute using the 2 (DIM_2D) or 3(DIM_3D) nearest states
///@todo compute using the projection on the nearest face

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  
  std::vector< Common::SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->getTrsList();

  const CFuint nbTrsInOtherMesh = trsList.size();

  CFreal minDistance = MathTools::MathConsts::CFrealMax();
  std::vector<CFuint> closestNodes(PhysicalModelStack::getActive()->getDim());
  for (CFuint iTrs = 0; iTrs < nbTrsInOtherMesh; ++iTrs) {

    Common::SafePtr<TopologicalRegionSet> trs = trsList[iTrs];

    const std::string trsName = trs->getName();
    if((trsName != "InnerCells") && (trsName != "InnerCells") && (trsName != "InnerFaces")){

      const CFuint nbOtherStates = trs->getNodesInTrs()->size();
      for(CFuint iOtherState=0; iOtherState < nbOtherStates; iOtherState++)
      {
        const CFuint otherNodeID = (*(trs->getNodesInTrs()))[iOtherState];
        _tempVector = *(otherNodes[otherNodeID]) - _coord;
        const CFreal distance = _tempVector.norm2();

        if(distance < minDistance)
        {

          if(closestNodes.size() == DIM_3D) closestNodes[2] = closestNodes[1];
          closestNodes[1] = closestNodes[0];
          closestNodes[0] = otherNodeID;
          minDistance = distance;
        }
      }
    }
  }

  //Compute the interpolated state from the 2 or 3 closest states
  CFreal sumOverDistance = 0;
  bool exact(false);
  for(CFuint iNState=0;iNState < closestNodes.size();iNState++)
  {
    const CFuint otherNodeID = closestNodes[iNState];
    _tempVector = *(otherNodes[otherNodeID]) - _coord;
    const CFreal distance = _tempVector.norm2();

    if(distance < MathTools::MathConsts::CFrealEps())
    {
      _interpolatedState = (otherStates[otherNodeID]);
      exact = true;
    }
    else{
      _interpolatedState += (otherStates[otherNodeID]/distance);
      sumOverDistance += 1./distance;
    }
  }
  if(!exact) _interpolatedState /= sumOverDistance;

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  LinearMeshInterpolatorFVMCC::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = LinearMeshInterpolator::needsSockets();

  result.push_back(&socket_otherNstates);
  result.push_back(&socket_otherGstates);
  result.push_back(&socket_otherNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
