#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"

#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/TrsNotFoundException.hh"

#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/VolumeBasedExtrapolator.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<VolumeBasedExtrapolator,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeModule>
VolumeBasedExtrapolatorProvider("VolumeBased");

//////////////////////////////////////////////////////////////////////////////

VolumeBasedExtrapolator::VolumeBasedExtrapolator(const std::string& name) :
  NodalStatesExtrapolator<CellCenterFVMData>(name),
  socket_volumes("volumes"),
  _weights(),
  _weightsStorage()
{
}

//////////////////////////////////////////////////////////////////////////////

VolumeBasedExtrapolator::~VolumeBasedExtrapolator()
{
}

//////////////////////////////////////////////////////////////////////////////

void VolumeBasedExtrapolator::extrapolateInAllNodes()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();
  cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D ||
   PhysicalModelStack::getActive()->getDim() == DIM_3D);

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    // reset to 0 the nodal states and weights
    nodalStates[iNode] = 0.;
    const CFuint nbNeighborStates = _neighborStates[iNode].size();

    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const State *const neighState = _neighborStates[iNode][iState];
      const CFreal weight = static_cast<CFreal>(_weights[iNode][iState]);
      nodalStates[iNode] += (*neighState)*weight;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolumeBasedExtrapolator::extrapolateInNodes
(const vector<Node*>& nodes)
{

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

  const CFuint nbNodes = nodes.size();
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const CFuint nodeID = nodes[iNode]->getLocalID();
    // reset to 0 the nodal states
    nodalStates[nodeID] = 0.;
    const CFuint nbNeighborStates = _neighborStates[nodeID].size();

    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const CFreal weight = static_cast<CFreal>(_weights[nodeID][iState]);
      nodalStates[nodeID] += (*_neighborStates[nodeID][iState])*weight;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolumeBasedExtrapolator::setup()
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  // first call the parent class
  NodalStatesExtrapolator<CellCenterFVMData>::setup();

  const CFuint nbNodes = nodes.size();
  cf_assert(PhysicalModelStack::getActive()->getDim() == 2 ||
   PhysicalModelStack::getActive()->getDim() == 3);

  // resize the storage of the weights
  CFuint weightsStorageSize = 0;
  for (CFuint i = 0; i < nbNodes; ++i) {
    weightsStorageSize += _neighborStates[i].size();
  }

  CFLog(VERBOSE, "Memory Usage before weights: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n" );

  // resize the table
  _weights.resize(nbNodes);
  _weightsStorage.resize(weightsStorageSize, 0.0);

  CFuint counter = 0;
  for (CFuint i = 0; i < nbNodes; ++i) {
    cf_assert(counter < weightsStorageSize);
    _weights[i]  = &_weightsStorage[counter];
    counter += _neighborStates[i].size();
  }

  CFLog(VERBOSE, "Memory Usage after weights: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n" );

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const CFuint nbNeighborStates = _neighborStates[iNode].size();

    CFreal sumWeight = 0.0;
    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      CFuint neighCellID;
      if(_neighborStates[iNode][iState]->isGhost())
         neighCellID = _ghost2InnerStates[_neighborStates[iNode][iState]->getLocalID()];
      else neighCellID = _neighborStates[iNode][iState]->getLocalID();

      const CFreal weight = 1./volumes[neighCellID];
      _weights[iNode][iState] = weight;
      sumWeight += weight;
    }

    // divide the weights by the sum of the weights
    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      _weights[iNode][iState] /= sumWeight;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolumeBasedExtrapolator::addBoundaryNeighbors
(vector< vector<CFuint> >& bFacesPerNode,
 CFMap<CFuint, CFint>& mapFaceTrsID)
{
  CFAUTOTRACE;

  DataHandle<CFint> trsID = socket_trsID.getDataHandle();
CFout <<"Here...\n";
  const CFuint nbGhostStates = socket_gstates.getDataHandle().size();
  _ghost2InnerStates.resize(nbGhostStates);

  // reset all the trsIDs to -1
  trsID = -1;

  vector< SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  vector<Common::SafePtr<TopologicalRegionSet> > trs;
  if (_trsPriorityList.size() > 0) {
    for (CFuint i = 0; i < _trsPriorityList.size(); ++i) {
      const std::string name = _trsPriorityList[i];
      bool nameFound = false;
      for (CFuint j = 0; j < trsList.size(); ++j) {
  if (name == trsList[j]->getName()) {
    trs.push_back(trsList[j]);
    nameFound = true;
    break;
  }
      }
      if (!nameFound) {
  throw TrsNotFoundException (FromHere(), name + " not found!!");
      }
    }
  }
  else {
    trs = trsList;
  }

  GeometricEntityPool<FaceTrsGeoBuilder> geoBuilder;
  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  geoBuilder.setup();

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  const CFint nbTRSs = trs.size();
  for (CFint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

    if (currTrs->getName() != "InnerCells" &&
        currTrs->getName() != "InnerFaces" &&
        currTrs->getName() != "PartitionFaces") {

      geoData.trs = currTrs;
      geoData.isBFace = true;

      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        CFLogDebugMed( "iFace = " << iFace << "\n");

        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity *const currFace = geoBuilder.buildGE();

        State* const gstate = currFace->getState(1);
        cf_assert(gstate != CFNULL);

        const vector<Node*>* const bnodes = currFace->getNodes();
        const CFuint nbNodesInFace = bnodes->size();
        for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
          const CFuint nodeID = (*bnodes)[iNode]->getLocalID();
          // set this face in the list of boundary faces referencing
          // this node
          bFacesPerNode[nodeID].push_back(currFace->getID());

          // if the trsID is different from -1 or the current one
          // nothing will happen : this is to prevent to distribute
          // contributions to CORNER boundary nodes from different TRSs
          if (trsID[nodeID] == -1) {
            // if the trsID hasn't been set yet, it is done and
            // the corresponding nodal state and sum of the weights are
            // reset to 0.
            trsID[nodeID] = iTRS;

      // register the name of the TRS with the corresponding ID
      _mapTrsNameToID.insert(currTrs->getName(), iTRS);
          }

          if (trsID[nodeID] == iTRS) {
            // this neighbor must be considered
            _neighborStates[nodeID].push_back(gstate);
            _ghost2InnerStates[gstate->getLocalID()] =  currFace->getState(0)->getLocalID();
          }

          mapFaceTrsID.insert(currFace->getID(), iTRS);
        }

        geoBuilder.releaseGE();
      }
    }
  }

  mapFaceTrsID.sortKeys();
  _mapTrsNameToID.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
