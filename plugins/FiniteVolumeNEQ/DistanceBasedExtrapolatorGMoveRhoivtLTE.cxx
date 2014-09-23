#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveRhoivtLTE.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveRhoivtLTE,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveRhoivtLTEProvider("DistanceBasedGMoveRhoivtLTE");

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtLTE::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<std::string> >
    ("TrsWithNoLTE", "list of the names of the TRS on which LTE has not to be applied");
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveRhoivtLTE::DistanceBasedExtrapolatorGMoveRhoivtLTE
(const std::string& name) :
  DistanceBasedExtrapolatorGMoveRhoivt(name),
  _xs(),
  _doNotApplyLTE()
{
  addConfigOptionsTo(this);
  
  _trsWithNoLTE = vector<std::string>();
  setParameter("TrsWithNoLTE",&_trsWithNoLTE);
}
      
//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveRhoivtLTE::~DistanceBasedExtrapolatorGMoveRhoivtLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtLTE::extrapolate()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  // reset to 0 the nodal states
  nodalStates[_currNodeID] = 0.;
  
  if (!_isNodeToPrescribe[_currNodeID]) { 
    applyInner(); 
  }
  else {
    (_doNotApplyLTE[_currNodeID]) ? 
      DistanceBasedExtrapolatorGMoveRhoivt::applyBC() : applyBC();
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtLTE::setup()
{
  DistanceBasedExtrapolatorGMoveRhoivt::setup();

  const CFuint nbSpecies = _library->getNbSpecies();
  cf_assert(nbSpecies > 0);
  _xs.resize(nbSpecies);
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  _doNotApplyLTE.resize(nodes.size());
  _doNotApplyLTE.assign(_doNotApplyLTE.size(), false);
  
  for (CFuint iTRS = 0; iTRS < _trsWithNoLTE.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_trsWithNoLTE[iTRS]);
    SafePtr<vector<CFuint> > nodesinTRS = trs->getNodesInTrs();
    const CFuint nbNodesInTRS = nodesinTRS->size();
    for (CFuint i = 0; i < nbNodesInTRS; ++i) {
      _doNotApplyLTE[(*nodesinTRS)[i]] = true;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtLTE::transformBack
(RealVector& nstate)
{
  DataHandle<CFint> trsID = socket_trsID.getDataHandle();

  // transform [y_i v p] into [rho_i v T Tv_i] variables, T and Tv_i
  // are already known because they are strongly imposed
  const CFuint nbSpecies = _library->getNbSpecies();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = nbSpecies + dim;

  // T can be imposed or coming from the radiative equilibrium condition
  CFreal T = 0.0;
  if ((!_radEquilibrium || _doNotApplyRadEq[_currNodeID]) && !runAdiabatic()) {
    if (!_nodalValuesIDFlags[TID]) {
	T = _wValues[dim];
    }
    else {
      T = getNodalValue(_orderedTrsList[trsID[_currNodeID]],TID,_currNodeID);
      
      // other vibrational temperatures (careful not to override Te)
      const CFuint nbTs = _library->getNbTempVib() - _library->getNbTe() + 1;
      for (CFuint i = 0; i < nbTs; ++i) {
	nstate[TID + i] = T;
      }
    }
  }
  else {
    T = nstate[TID];
  }
    
  CFreal p = _nodalPressure; // p is the reconstructed pressure
  
  const RealVector& refData = _model->getReferencePhysicalData();
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*refData[EulerTerm::T];
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = nstate[ie]; // reconstructed mass fractions
  }

  // the elemental composition Xe corresponding to the extrapolated
  // species composition is computed
  _library->setElementXFromSpeciesY(_ys);
  
  // the LTE composition corresponding to (T,p,Xe)
  _library->setComposition(Tdim,pdim,&_xs);
  
  // the density corresponding to the LTE composition
  CFreal* TeDim = (_library->getNbTempVib() > 0) ? &nstate[TID+1] : CFNULL;
  const CFreal rhoEq = _library->density(Tdim, pdim, TeDim)/refData[EulerTerm::RHO];
  
  // the mass fractions corresponding to the equilibrium conditions
  _library->getSpeciesMassFractions(_xs,_ys);
  
  // rho_i are consistent with the wall temperature, the reconstructed pressure,
  // the assumption of LTE  composition at the wall
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    // partial densities must be positive
    nstate[ie] = max(0.0, rhoEq*_ys[ie]);
  }

  // velocity components and temperature are simply overridden
  for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
    if (_wValuesIdx[i] < TID)  {
      nstate[_wValuesIdx[i]] = _wValues[i];
    }
    else {
      // temperatures are assigned to computed temperature
      nstate[_wValuesIdx[i]] = T;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtLTE::applyBC()
{    
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  const Node& currNode = *nodes[_currNodeID];
  const CFuint nbNeighborStates = _neighborStates[_currNodeID].size();
  const CFuint nbEqs =  PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbSpecies = _library->getNbSpecies();
   
  // update the weights that might have changed because of the ghost state repositioning
  updateWeights(_currNodeID, currNode);
  
  _nodalPressure = 0.0;
  CFreal sumDr = 0.0;
  
  for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
    const CFreal weight = static_cast<CFreal>(_weights[_currNodeID][iState]);
    const State *const neighState = _neighborStates[_currNodeID][iState];
    
    transform(*neighState, _tmpState);
    
    // only internal States are considered
    // in fact I'm here only interested in extrapolating p, y_i, Te
    CFreal invR = 0.0;
    if (!neighState->isGhost()) {
      invR = 1./MathFunctions::getDistance(neighState->getCoordinates(),
					   *nodes[_currNodeID]);
      sumDr += invR;
    }
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      if (iVar < nbSpecies) {
	// in this case the weights have to be locally recomputed from
	// only internal states
	if (!neighState->isGhost()) {
	  // y_i are extrapolated only from inside the domain
	  nodalStates[_currNodeID][iVar] += _tmpState[iVar]*invR;
	}
      }
      else {
	// other variables are extrapolated from ghost and inner states
	nodalStates[_currNodeID][iVar] += _tmpState[iVar]*weight;
      }
    }
    _nodalPressure += _pressure*weight;
  }
  
  for (CFuint iVar = 0; iVar < nbSpecies; ++iVar) {
    nodalStates[_currNodeID][iVar] /= sumDr;
  }
  
  transformBack(nodalStates[_currNodeID]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
