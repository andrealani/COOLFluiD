#include "FiniteVolume/DistanceBasedExtrapolatorGMoveCoupled.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/State.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveCoupled,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeModule>
DistanceBasedExtrapolatorGMoveCoupledProvider("DistanceBasedGMoveCoupled");

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCoupled::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("CoupledValuesIdx","Indices of the variables transfered through coupling");

  options.addConfigOption< vector<CFreal> >
    ("DefaultCoupledValues","Default values for the variables transfered through coupling");

  options.addConfigOption< vector<std::string> >
    ("Interfaces","Names of the Interfaces for the transfer");

  options.addConfigOption< CFuint >
    ("DefaultIterations","Default number of iterations to use the default values");

}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveCoupled::DistanceBasedExtrapolatorGMoveCoupled
(const std::string& name) :
  DistanceBasedExtrapolatorGMove(name),
  _sockets(),
  _isSetIndex(false)
{

  addConfigOptionsTo(this);

  _wCoupledValuesIdx = vector<CFuint>();
  setParameter("CoupledValuesIdx",&_wCoupledValuesIdx);

  _defaultCoupledValues = vector<CFreal>();
  setParameter("DefaultCoupledValues",&_defaultCoupledValues);

  _interfaceNames = vector<std::string>();
  setParameter("Interfaces",&_interfaceNames);

  _defaultIterations = 10;
  setParameter("DefaultIterations",&_defaultIterations);

}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveCoupled::~DistanceBasedExtrapolatorGMoveCoupled()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCoupled::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DistanceBasedExtrapolatorGMove::configure(args);

  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  
  cf_assert(currentSubSystem == SubSystemStatusStack::getCurrentName());
  cf_assert(_interfaceNames.size() == _trsName.size());
  cf_assert(_trsName.size() > 0);

  for(CFuint iTRS = 0; iTRS < _trsName.size(); iTRS++)
  {
    const std::string trsName = _trsName[iTRS];
    const std::string interfaceName = _interfaceNames[iTRS];

    const std::string baseSocketName =
      "COUPLING_" + interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";
    std::string socketName = baseSocketName + "ISACCEPTED";
    _sockets.createSocketSink<CFreal>(socketName);
    socketName = baseSocketName + "DATA";
    _sockets.createSocketSink<RealVector>(socketName);
  }
}


//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCoupled::setup()
{
  CFAUTOTRACE;

  cf_assert(_trsName.size() > 0);

  DistanceBasedExtrapolatorGMove::setup();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();

  //Create a Map to get back the index of the node in the TRS list from its LocalID
  _trsNodeIDMap.resize(_trsName.size());
  for (CFuint iName = 0; iName < _trsName.size(); ++iName) {
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

      if (currTrs->getName() == _trsName[iName]) {

        Common::SafePtr< vector<CFuint> > const trsNodes = currTrs->getNodesInTrs();
        const CFuint nbNodesInTRS = trsNodes->size();
        for (CFuint iNode = 0; iNode < nbNodesInTRS; ++iNode) {
          const CFuint nodeID = (*trsNodes)[iNode];
          _trsNodeIDMap[iName].insert(nodeID, iNode);
        }
      }
    }
    _trsNodeIDMap[iName].sortKeys();
  }

  //check that prescribed values and coupled
  //values idx are different!!!
  cf_assert(_wCoupledValuesIdx.size() > 0);
  cf_assert(_wValuesIdx.size() > 0);
  for(CFuint i=0;i<_wCoupledValuesIdx.size();++i){
    for(CFuint j=0;j < _wValuesIdx.size();++j){
      cf_assert(_wCoupledValuesIdx[i] != _wValuesIdx[j]);
    }
  }

  //Resize the coupled values vector
  _coupledValues.resize(_wCoupledValuesIdx.size());

  cf_assert(_defaultCoupledValues.size() == _wCoupledValuesIdx.size());

}


//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCoupled::setIndex()
{
  CFAUTOTRACE;

  cf_assert(_trsName.size() > 0);

  _coupledDataID.resize(_trsName.size());
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();

  for (CFuint iName = 0; iName < _trsName.size(); ++iName) {
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

      if (currTrs->getName() == _trsName[iName]) {
        Common::SafePtr< vector<CFuint> > const trsNodes = currTrs->getNodesInTrs();
        const CFuint nbNodesInTRS = trsNodes->size();

        const std::string trsName = currTrs->getName();
        const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
        const std::string nameSpace = getMethodData().getNamespace();
        const std::string interfaceName = _interfaceNames[iName];
        const std::string baseSocketName =
          "COUPLING_" + interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";
        std::string socketName = baseSocketName + "ISACCEPTED";
        DataHandle<CFreal> isAccepted =
          _sockets.getSocketSink<CFreal>(socketName)->getDataHandle();

        cf_assert(isAccepted.size() == nbNodesInTRS);
        _coupledDataID[iName].resize(nbNodesInTRS);
        CFuint idx = 0;
        for(CFuint iNode=0; iNode < nbNodesInTRS; ++iNode)
        {
          if(isAccepted[iNode] >= 0.)
          {
            _coupledDataID[iName][iNode] = idx;
            idx++;
          }
          else{
            _coupledDataID[iName][iNode] = -1;
          }
        }
      }
    }
  }
  _isSetIndex = true;
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCoupled::extrapolateInAllNodes()
{
  CFAUTOTRACE;

  if(!_isSetIndex) setIndex();

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();
  cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D ||
	 PhysicalModelStack::getActive()->getDim() == DIM_3D);

// unused //  const CFuint nbEqs =  PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const CFuint nbNeighborStates = _neighborStates[iNode].size();

    // reset to 0 the nodal states
    nodalStates[iNode] = 0.;
    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const CFreal weight = static_cast<CFreal>(_weights[iNode][iState]);
      const State *const neighState = _neighborStates[iNode][iState];
      nodalStates[iNode] += (*neighState)*weight;
    }
  }

  //For the boundary nodes
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();

// unused //  CFuint count = 0;
  for (CFuint iName = 0; iName < _trsName.size(); ++iName) {
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

      if (currTrs->getName() == _trsName[iName]) {
        Common::SafePtr< vector<CFuint> > const trsNodes = currTrs->getNodesInTrs();
        const CFuint nbTrsNodes = trsNodes->size();
        const std::string trsName = currTrs->getName();
        const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
        const std::string nameSpace = getMethodData().getNamespace();
        const std::string interfaceName = _interfaceNames[iName];
        const std::string baseSocketName =
          "COUPLING_" + interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";
        std::string socketName = baseSocketName + "ISACCEPTED";
        DataHandle<CFreal> isAccepted =
          _sockets.getSocketSink<CFreal>(socketName)->getDataHandle();

        socketName = baseSocketName + "DATA";
        DataHandle<RealVector> interfaceData =
          _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();
        cf_assert(isAccepted.size() == nbTrsNodes);
        for(CFuint iNode=0; iNode < nbTrsNodes; ++iNode)
        {
          const CFuint nodeID = nodes[(*trsNodes)[iNode]]->getLocalID();
          if (_isNodeToPrescribe[nodeID]) {
            for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
              nodalStates[nodeID][_wValuesIdx[i]] = _wValues[i];
            }

            //Get the value to impose at node iNode
            if((isAccepted[iNode]>=0.) && SubSystemStatusStack::getActive()->getNbIter() > _defaultIterations){
              _coupledValues = interfaceData[_coupledDataID[iName][iNode]];
              for (CFuint i = 0; i < _wCoupledValuesIdx.size(); ++i) {
                nodalStates[nodeID][_wCoupledValuesIdx[i]] = _coupledValues[i];
              }
            }
            else{
              for (CFuint i = 0; i < _wCoupledValuesIdx.size(); ++i) {
                nodalStates[nodeID][_wCoupledValuesIdx[i]] = _defaultCoupledValues[i];
              }
            }
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCoupled::extrapolateInNodes
(const vector<Node*>& nodes)
{
  cf_assert(false);
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
