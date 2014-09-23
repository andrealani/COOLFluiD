#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveRhoivtCat.hh"
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

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveRhoivtCat,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveRhoivtCatProvider("DistanceBasedGMoveRhoivtCat");

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtCat::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<std::string> >
    ("TrsWithNoCat", "list of the names of the TRS on which Cat has not to be applied");
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveRhoivtCat::DistanceBasedExtrapolatorGMoveRhoivtCat
(const std::string& name) :
  DistanceBasedExtrapolatorGMoveRhoivt(name),
  _doNotApplyCat()
{
  addConfigOptionsTo(this);
  
  _trsWithNoCat = vector<std::string>();
  setParameter("TrsWithNoCat",&_trsWithNoCat);
}
      
//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveRhoivtCat::~DistanceBasedExtrapolatorGMoveRhoivtCat()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtCat::extrapolate()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  // reset to 0 the nodal states
  nodalStates[_currNodeID] = 0.;
  
  if (!_isNodeToPrescribe[_currNodeID]) { 
    applyInner(); 
  }
  else {
    (_doNotApplyCat[_currNodeID]) ? 
      DistanceBasedExtrapolatorGMoveRhoivt::applyBC() : applyBC();
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtCat::setup()
{
  DistanceBasedExtrapolatorGMoveRhoivt::setup();
  
  cf_assert( _library->getNbSpecies() > 0);
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  _doNotApplyCat.resize(nodes.size());
  _doNotApplyCat.assign(_doNotApplyCat.size(), false);
  
  for (CFuint iTRS = 0; iTRS < _trsWithNoCat.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_trsWithNoCat[iTRS]);
    SafePtr<vector<CFuint> > nodesinTRS = trs->getNodesInTrs();
    const CFuint nbNodesInTRS = nodesinTRS->size();
    for (CFuint i = 0; i < nbNodesInTRS; ++i) {
      _doNotApplyCat[(*nodesinTRS)[i]] = true;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivtCat::transformBack(RealVector& nstate)
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
    _ys[ie] = nstate[ie]; 
    // reconstructed mass fractions: this is inconsistent with BC values
    // if ghost states are moved from original positions, but this should not
    // be the case if we start from a developed solution
  }

  // density on the wall given from equation of state
  _library->setSpeciesFractions(_ys);
  /*const CFreal coeff = _library->getRgas()/_library->getMMass();
  const CFreal rhoWall = pdim/(Tdim*coeff)/refData[EulerTerm::RHO];
  
  // rho_i are consistent with the wall temperature, the reconstructed pressure,
  // the assumption of Cat  composition at the wall
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    // partial densities must be positive
    nstate[ie] = max(0.0, rhoWall*_ys[ie]);
  }*/

 
  // rho_i are consistent with the wall temperature and the reconstructed pressure
  CFreal* TeDim = (_library->getNbTempVib() > 0) ? &nstate[TID+1] : CFNULL;
  const CFreal rho = _library->density(Tdim, pdim, TeDim)/refData[EulerTerm::RHO];
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    nstate[ie] = max(0.0, rho*_ys[ie]);
  }

  // velocity components and temperature are simply overridden
  for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
    nstate[_wValuesIdx[i]] = _wValues[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
