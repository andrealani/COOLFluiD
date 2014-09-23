#include "DistanceBasedExtrapolatorGMove.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Common/BadValueException.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

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

MethodStrategyProvider<DistanceBasedExtrapolatorGMove,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeModule>
distanceBasedExtrapolatorGMoveProvider("DistanceBasedGMove");

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMove::DistanceBasedExtrapolatorGMove
(const std::string& name) :
  DistanceBasedExtrapolator<CellCenterFVMData>(name),
  _isNodeToPrescribe(),
  _doNotApplyRadEq()
{
  addConfigOptionsTo(this);
  
  _radTID = 10000000;
  setParameter("RadTID",&_radTID);
  
  _radEquilibrium = false;
  setParameter("RadEquilibrium",&_radEquilibrium);

  _wValuesIdx = vector<CFuint>();
  setParameter("ValuesIdx",&_wValuesIdx);

  _wValues = vector<CFreal>();
  setParameter("Values",&_wValues);
  
  _trsWithNoRadEq = vector<std::string>();
  setParameter("TrsWithNoRadEquilibrium",&_trsWithNoRadEq);
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMove::~DistanceBasedExtrapolatorGMove()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMove::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint > ("RadTID", "Temperature ID for controlling radiative equilibrium");
  options.addConfigOption< bool, Config::DynamicOption<> >
    ("RadEquilibrium", "Flag to activate radiative equilibrium");
  options.addConfigOption< vector<CFuint> >
    ("ValuesIdx","Indices of the prescribed values");
  options.addConfigOption< vector<CFreal> >
    ("Values","Prescribed values"); 
  options.addConfigOption< vector<std::string> > 
    ("TrsWithNoRadEquilibrium", 
     "list of the names of TRSs on which Radiative Equilibrium has not to be applied");
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMove::setup()
{
  DistanceBasedExtrapolator<CellCenterFVMData>::setup();

  cf_assert(_trsName.size() > 0);
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFint> trsID = socket_trsID.getDataHandle();

  _isNodeToPrescribe.resize(nodes.size());
  _isNodeToPrescribe = false;

  const CFuint nbTRSs = _orderedTrsList.size();

  CFuint count = 0;
  for (CFuint iName = 0; iName < _trsName.size(); ++iName) {
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<TopologicalRegionSet> currTrs = _orderedTrsList[iTRS];

      if (currTrs->getName() == _trsName[iName]) {
	Common::SafePtr<std::vector<CFuint> > trsNodes = currTrs->getNodesInTrs();
	const CFuint nbTrsNodes = trsNodes->size();

        for (CFuint i = 0; i < nbTrsNodes; ++i) {
          const CFuint nodeID = (*trsNodes)[i];

          // only nodes which are flagged to belong to this TRS can be prescribed
          if (trsID[nodeID] == (CFint) iTRS) {
	    _isNodeToPrescribe[nodeID] = true;
          }
	}
	count++;
	break;
      }
    }
  }

  checkTRSList(count);
  
  _doNotApplyRadEq.resize(nodes.size());
  _doNotApplyRadEq.assign(_doNotApplyRadEq.size(), false);
  
  for (CFuint iTRS = 0; iTRS < _trsWithNoRadEq.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_trsWithNoRadEq[iTRS]);
    SafePtr<vector<CFuint> > nodesinTRS = trs->getNodesInTrs();
    const CFuint nbNodesInTRS = nodesinTRS->size();
    for (CFuint i = 0; i < nbNodesInTRS; ++i) {
      _doNotApplyRadEq[(*nodesinTRS)[i]] = true;
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMove::checkTRSList(const CFuint count)
{
  if (count != _trsName.size()) {
    throw BadValueException
      (FromHere(),"DistanceBasedExtrapolatorGMove::setup() => TRS name not found");
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMove::applyBC()
{
  if (_radTID == 10000000 && _radEquilibrium) {
    CFLog(WARN, "DistanceBasedExtrapolatorGMove::setup() => RadTID not correctly set !!!");
    cf_assert(_radTID < 10000000 && _radEquilibrium);
  }
  
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle<CFint> trsID = socket_trsID.getDataHandle();

  for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
    const CFuint idx = _wValuesIdx[i];
    if (idx < _radTID) {
      nodalStates[_currNodeID][idx] = (!_nodalValuesIDFlags[idx]) ?
	_wValues[i] : getNodalValue(_orderedTrsList[trsID[_currNodeID]],idx,_currNodeID);
    }
    else {
      // temperature is fixed only if radiative equilibrium is not in use
      if ((!_radEquilibrium || _doNotApplyRadEq[_currNodeID]) && !runAdiabatic()) {
	nodalStates[_currNodeID][idx] = (!_nodalValuesIDFlags[idx]) ?
	  _wValues[i] : getNodalValue(_orderedTrsList[trsID[_currNodeID]],idx,_currNodeID);
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMove::extrapolate()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

  // reset to 0 the nodal states
  nodalStates[_currNodeID] = 0.;
  applyInner();

  if (_isNodeToPrescribe[_currNodeID]) {
    applyBC();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
