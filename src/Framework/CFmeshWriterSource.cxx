// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/ConnectivityTable.hh"
#include "Common/CFLog.hh"

#include "Framework/MeshData.hh"
#include "Framework/CFmeshWriterSource.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CFmeshWriterSource::CFmeshWriterSource() :
  _nodes(CFNULL),
  _states(CFNULL),
  _cellNodes(CFNULL),
  _cellStates(CFNULL),
  _nodalExtraDataSinkSockets(),
  _stateExtraDataSinkSockets()
{
}

//////////////////////////////////////////////////////////////////////////////

CFmeshWriterSource::~CFmeshWriterSource()
{
  releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::releaseMemory()
{
  BaseCFMeshFileSource::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::setMeshData()
{
  // set the handles for the needed storages
  _nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  _states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  _cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  _cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  
  // set those other data needed during the writing
  copyMeshData();
  copyTrsData();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::setExtraDataSockets(Common::SafePtr<DynamicDataSocketSet<> > dynamicSocket)
{
  // set the sockets for the needed storages
  dynamicSockets = dynamicSocket;

  copyExtraData();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::copyExtraData()
{
  _nbExtraNodalVars = _extraNVarMap->size();
  _nbExtraStateVars = _extraSVarMap->size();
  _nbExtraVars = _extraVarMap->size();

  _extraNodalVarNames.resize(_nbExtraNodalVars);
  _extraStateVarNames.resize(_nbExtraStateVars);
  _extraVarNames.resize(_nbExtraVars);

  _extraNodalVarStrides.resize(_nbExtraNodalVars);
  _extraStateVarStrides.resize(_nbExtraStateVars);
  _extraVarStrides.resize(_nbExtraVars);

  CFuint totalNbExtraNodalVars = 0;
  for(CFuint iVar=0; iVar < _nbExtraNodalVars; iVar++)
  {
    // this is dangerous because the order at this moment
    // can be sorted or not ...
    _extraNodalVarNames[iVar]   = (*_extraNVarMap)[iVar].first;
    _extraNodalVarStrides[iVar] = (*_extraNVarMap)[iVar].second;
    totalNbExtraNodalVars += _extraNodalVarStrides[iVar];
  }

  CFuint totalNbExtraStateVars = 0;
  for(CFuint iVar=0; iVar < _nbExtraStateVars; iVar++)
  {
    // this is dangerous because the order at this moment
    // can be sorted or not ...
    _extraStateVarNames[iVar]   = (*_extraSVarMap)[iVar].first;
    _extraStateVarStrides[iVar] = (*_extraSVarMap)[iVar].second;
    totalNbExtraStateVars += _extraStateVarStrides[iVar];
  }
  
  CFuint totalNbExtraVars = 0;
  for(CFuint iVar=0; iVar < _nbExtraVars; iVar++)
  {
    // this is dangerous because the order at this moment
    // can be sorted or not ...
    _extraVarNames[iVar]   = (*_extraVarMap)[iVar].first;
    _extraVarStrides[iVar] = (*_extraVarMap)[iVar].second;
    totalNbExtraVars += _extraVarStrides[iVar];
  }

  _extraStateVector.resize(totalNbExtraStateVars);
  _extraNodeVector.resize(totalNbExtraNodalVars);
  _extraVector.resize(totalNbExtraVars);

}

//////////////////////////////////////////////////////////////////////////////

const RealVector* CFmeshWriterSource::getPastNode(const CFuint nodeID) const
{
  return (dynamicSockets->getSocketSink<Node*>("pastNodes")->getDataHandle())[nodeID];
}

//////////////////////////////////////////////////////////////////////////////

const RealVector* CFmeshWriterSource::getPastState(const CFuint stateID) const
{
  return (dynamicSockets->getSocketSink<State*>("pastStates")->getDataHandle())[stateID];
}

//////////////////////////////////////////////////////////////////////////////

const RealVector* CFmeshWriterSource::getInterNode(const CFuint nodeID) const
{
  return (dynamicSockets->getSocketSink<Node*>("interNodes")->getDataHandle())[nodeID];
}

//////////////////////////////////////////////////////////////////////////////

const RealVector* CFmeshWriterSource::getInterState(const CFuint stateID) const
{
  return (dynamicSockets->getSocketSink<State*>("interStates")->getDataHandle())[stateID];
}
//////////////////////////////////////////////////////////////////////////////

RealVector& CFmeshWriterSource::getExtraStateValues(const CFuint stateID)
{
  SafePtr<vector<CFuint> > strides = getExtraStateVarStrides();
  CFuint count = 0;

  cf_assert(_stateExtraDataSinkSockets.size() == _nbExtraStateVars);
  for(CFuint iVar=0; iVar < _nbExtraStateVars; ++iVar)
  {
    DataHandle<CFreal> extraVars = _stateExtraDataSinkSockets[iVar]->getDataHandle();
    const CFuint currStride = (*strides)[iVar];
    const CFuint start = stateID*currStride;

    cf_assert(start < extraVars.size());
    for (CFuint i = 0; i < currStride; ++i, ++count) {
      _extraStateVector[count] = extraVars[start+i];
    }
  }

  return _extraStateVector;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& CFmeshWriterSource::getExtraNodalValues(const CFuint nodeID)
{
  SafePtr<vector<CFuint> > strides = getExtraNodalVarStrides();
  CFuint count = 0;

  for(CFuint iVar=0; iVar < _nbExtraNodalVars; ++iVar)
  {
    DataHandle<CFreal> extraVars = _nodalExtraDataSinkSockets[iVar]->getDataHandle();
    const CFuint currStride = (*strides)[iVar];
    const CFuint start = nodeID*currStride;
    for (CFuint i = 0; i < currStride; ++i, ++count) {
      _extraNodeVector[count] = extraVars[start+i];
    }
  }

  return _extraNodeVector;
}
    
//////////////////////////////////////////////////////////////////////////////

RealVector& CFmeshWriterSource::getExtraValues()
{
  SafePtr<vector<CFuint> > strides = getExtraVarStrides();
  CFuint count = 0;
  
  for(CFuint iVar=0; iVar < _nbExtraVars; ++iVar)
  {
    DataHandle<CFreal> extraVars = _extraDataSinkSockets[iVar]->getDataHandle();
    const CFuint currStride = (*strides)[iVar];
    const CFuint start = 0;
    for (CFuint i = 0; i < currStride; ++i, ++count) {
      _extraVector[count] = extraVars[start+i];
    }
  }
  
  return _extraVector;
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::copyMeshData()
{
  _isWithSolution = true;

  _dimension = PhysicalModelStack::getActive()->getDim();

  _nbEquations = PhysicalModelStack::getActive()->getNbEq();

  /// @todo fix this for parallel mesh reading
  _nbUpdatableNodes = _nodes.size();
  _nbNonUpdatableNodes = 0;
  _nbUpdatableStates = _states.size();
  _nbNonUpdatableStates = 0;

  cf_assert(_cellStates->nbRows() == _cellNodes->nbRows());
  _nbElements = _cellStates->nbRows();

  SafePtr<vector<ElementTypeData> > elementTypeData =
    MeshDataStack::getActive()->getElementTypeData();

  // copy the element type data from the ones stored in MeshData
  _elementTypeData = *elementTypeData;

  _nbElementTypes = _elementTypeData.size();

  _geometricPolyOrder = static_cast<CFPolyOrder::Type>
    (_elementTypeData[0].getGeoOrder());

  _solutionPolyOrder = static_cast<CFPolyOrder::Type>
    (_elementTypeData[0].getSolOrder());
}
//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::copyTrsData()
{
  CFAUTOTRACE;

  vector< Common::SafePtr<TopologicalRegionSet> > alltrs = MeshDataStack::getActive()->getTrsList();

  CFuint countTRS = 0;
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = alltrs.begin(); itrs != alltrs.end(); ++itrs) {
    if (!(*itrs)->hasTag("writable")) ++countTRS;
  }
    
  const CFuint nbTRSs = alltrs.size() - countTRS;
  _nbTRSs = nbTRSs;
  _geoConn.resize(nbTRSs);
  _nbTRs.resize(nbTRSs);
  
  // only boundary TRS are listed
  CFuint iTRS = 0;
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator itr;
  for (itr = alltrs.begin(); itr != alltrs.end(); ++itr) {
    SafePtr<TopologicalRegionSet> trs = *itr;
    const std::string nameTRS = trs->getName();

    if (trs->hasTag("writable")) {
      _nameTRS.push_back(nameTRS);
      const CFuint nbTRsInTRS = trs->getNbTRs();
      _geoConn[iTRS].resize(nbTRsInTRS);
      _nbTRs[iTRS] = nbTRsInTRS;

      vector<CFuint> nbGeosInTR(0);
      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
        const CFuint nbGeos = (**itr)[iTR]->getLocalNbGeoEnts();
        nbGeosInTR.push_back(nbGeos);
        _geoConn[iTRS][iTR].resize(nbGeos);
      }

      _nbGeomEntsPerTR.push_back(nbGeosInTR);
      _geomType.push_back(CFGeoEnt::FACE);

      // geometric entities are ordered sequentially, TR after TR, in the TRS
      CFuint geoID = 0;
      for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
        SafePtr<TopologicalRegion> tr = (*trs)[iTR];
        const CFuint nbGeos = tr->getLocalNbGeoEnts();

        for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo, ++geoID) {
          const CFuint nbGeoNodes = trs->getNbNodesInGeo(geoID);
          CFuint nbGeoStates = trs->getNbStatesInGeo(geoID);
	  
          /// @todo gory fix (to see if this is cell centered FVM check if the   solution polyorder is 0)
          // if we are dealing with FVM cell centered meshes, boundary faces
          // have only one VALID state (the second one is a ghost-state)
          // and cannot be taken into account while outputting the mesh
          if (_solutionPolyOrder == 0 && nbGeoStates > 0) {
	    nbGeoStates = 1;
	  }
	  
          _geoConn[iTRS][iTR][iGeo].first.resize(nbGeoNodes);
          _geoConn[iTRS][iTR][iGeo].second.resize(nbGeoStates);

          for (CFuint iNode = 0; iNode < nbGeoNodes; ++iNode) {
            _geoConn[iTRS][iTR][iGeo].first[iNode] = trs->getNodeID(geoID,iNode);
          }
          for (CFuint iState = 0; iState < nbGeoStates; ++iState) {
            _geoConn[iTRS][iTR][iGeo].second[iState] = trs->getStateID(geoID,iState);
          }
        }
      }

      // only if the current TRS is writable increment the counter
      ++iTRS;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::prepareNodalExtraVars()
{
  // here the order of extra vars is the order imposed by the user in the
  // input CFcase file
  _nodalExtraDataSinkSockets.resize(_nbExtraNodalVars);
  for(CFuint iVar = 0; iVar < _nbExtraNodalVars;iVar++) {
    std::string socketName = (*_extraNVarMap)[iVar].first;
    _nodalExtraDataSinkSockets[iVar] = dynamicSockets->getSocketSink<CFreal>(socketName);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::prepareStateExtraVars()
{
  // here the order of extra vars is the order imposed by the user in the
  // input CFcase file
  _stateExtraDataSinkSockets.resize(_nbExtraStateVars);
  for(CFuint iVar = 0; iVar < _nbExtraStateVars;iVar++) {
    std::string socketName = (*_extraSVarMap)[iVar].first;
    _stateExtraDataSinkSockets[iVar] = dynamicSockets->getSocketSink<CFreal>(socketName);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterSource::prepareExtraVars()
{
  // here the order of extra vars is the order imposed by the user in the
  // input CFcase file
  _extraDataSinkSockets.resize(_nbExtraVars);
  for(CFuint iVar = 0; iVar < _nbExtraVars;iVar++) {
    std::string socketName = (*_extraVarMap)[iVar].first;
    _extraDataSinkSockets[iVar] = dynamicSockets->getSocketSink<CFreal>(socketName);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
