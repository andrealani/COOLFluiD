#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveRhoivt.hh"

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

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveRhoivt,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveRhoivtProvider("DistanceBasedGMoveRhoivt");

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivt::defineConfigOptions
(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveRhoivt::DistanceBasedExtrapolatorGMoveRhoivt
(const std::string& name) :
  DistanceBasedExtrapolatorGMove(name),
  _library(CFNULL),
  _model(CFNULL),
  _tmpState(),
  _ys(),
  _mmasses(),
  _nodalPressure(0.),
  _pressure(0.)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveRhoivt::~DistanceBasedExtrapolatorGMoveRhoivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivt::extrapolate()
{
  CFLog(DEBUG_MED, "DistanceBasedExtrapolatorGMoveRhoivt::extrapolate()\n");

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  // reset to 0 the nodal states
  nodalStates[_currNodeID] = 0.;
  
  (!_isNodeToPrescribe[_currNodeID]) ? applyInner() : applyBC();
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivt::applyBC()
{
  CFLog(DEBUG_MED, "DistanceBasedExtrapolatorGMoveRhoivt::applyBC()\n");
  
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const CFuint nbNeighborStates = _neighborStates[_currNodeID].size();
  
  const Node& currNode = *nodes[_currNodeID];
  updateWeights(_currNodeID, currNode);
  
  _nodalPressure = 0.0;
  for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
    const CFreal weight = static_cast<CFreal>(_weights[_currNodeID][iState]);
    const State *const neighState = _neighborStates[_currNodeID][iState];
    
    transform(*neighState, _tmpState);
    nodalStates[_currNodeID] += _tmpState*weight;
    _nodalPressure += _pressure*weight;
  }
  transformBack(nodalStates[_currNodeID]);
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivt::setup()
{
  CFLog(DEBUG_MED, "DistanceBasedExtrapolatorGMoveRhoivt::setup()\n");
  
  DistanceBasedExtrapolatorGMove::setup();

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>(); 
  cf_assert(_library.isNotNull());

  _model = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm> >();
  cf_assert(_model.isNotNull());

  _tmpState.resize(PhysicalModelStack::getActive()->getNbEq());
  _ys.resize(_library->getNbSpecies());
  
  _mmasses.resize(_library->getNbSpecies());
  _library->getMolarMasses(_mmasses);
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivt::transform(const RealVector& in,
						     RealVector& out)
{
  CFLog(DEBUG_MED, "DistanceBasedExtrapolatorGMoveRhoivt::transform()\n");
  
  // transform [rho_i v T Tv_i] variables into [y_i v T Tv_i]
  const CFuint nbSpecies = _library->getNbSpecies();

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += in[ie];
  }

  // Set the species
  const CFreal ovRho = 1./rho;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    out[ie] = in[ie]*ovRho;
    _ys[ie] = out[ie];
  }

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!!
  _library->setSpeciesFractions(_ys);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = nbSpecies + dim;
  const RealVector& refData = _model->getReferencePhysicalData();
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal T = in[TID];  
  CFreal Tdim = T*refData[EulerTerm::T];
  CFuint tvID = TID + 1;
  const CFuint nbTv = _library->getNbTempVib();
  for  (CFuint i = 0; i < nbTv; ++i, ++tvID) {
    out[tvID] = in[tvID];
  }
  out[TID] = T;
  
  // AL: this needs to be fixed for multi-temperature with MUTATION++
  CFreal* rhoi = &const_cast<RealVector&>(in)[0];
  if (nbTv > 0) {
    CFreal* TTv = &const_cast<RealVector&>(in)[TID];
    _library->setState(rhoi, TTv);
  }
  else {
    cf_assert(nbTv == 0);
    _library->setState(rhoi, &Tdim);
  }
  
  CFreal* tVec = (nbTv == 0) ? CFNULL : &(const_cast<RealVector&>(in).ptr()[TID+1]);
  CFreal pdim = _library->pressure(rhodim, Tdim, tVec);
  for  (CFuint i = 0; i < dim; ++i) {
    const CFuint iVar = nbSpecies + i;
    out[iVar] = in[iVar];
  }

  _pressure = pdim/refData[EulerTerm::P];
  
  if (_library->getNbTe() == 1) {
    out[tvID] = in[tvID];
  } 
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveRhoivt::transformBack(RealVector& nstate)
{
  CFLog(DEBUG_MED, "DistanceBasedExtrapolatorGMoveRhoivt::transformBack()\n");
  
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
      T = _wValues[dim]; // here you are assuming that velocity components come first 
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
    cf_assert(T > 0.);	
  }
  
  CFreal p = _nodalPressure; // p is the reconstructed pressure

  const RealVector& refData = _model->getReferencePhysicalData();
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*refData[EulerTerm::T];
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = nstate[ie]; // reconstructed mass fractions
  }

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other
  // thermodynamic quantity !!!
  _library->setSpeciesFractions(_ys);
  
  // rho_i are consistent with the wall temperature and the reconstructed pressure
  CFreal* TeDim = (_library->getNbTempVib() > 0) ? &nstate[TID+1] : CFNULL;
  const CFreal rho = _library->density(Tdim, pdim, TeDim)/refData[EulerTerm::RHO];
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    nstate[ie] = max(0.0, rho*_ys[ie]);
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

void DistanceBasedExtrapolatorGMoveRhoivt::setSpeciesVariables(RealVector& nstate)
{
  CFLog(DEBUG_MED, "DistanceBasedExtrapolatorGMoveRhoivt::setSpeciesVariables()\n");
  
  // nstate is in [rho_i v T Tv] variables
  // we want to recalculate p_i from _nodalPressure keeping mass fractions _ys constant 
  
  // we want to keep the same y_i and recalculate rho in function of p 
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbSpecies = _library->getNbSpecies();
  
  CFreal ovRho = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    ovRho += nstate[i];
  }    
  ovRho = 1./ovRho;
  
  const CFreal RT = _library->getRgas()*nstate[nbSpecies + dim];
  // assumption here: second temperature is electron temperature
  const CFreal RTe = (_library->presenceElectron()) ? _library->getRgas()*nstate[nbSpecies + dim + 1] : RT;
  CFreal RTovM = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _ys[i] = nstate[i]*ovRho;
    
    (i != 0) ? 
      RTovM += _ys[i]/_mmasses[i]*RT :
      RTovM += _ys[0]/_mmasses[0]*RTe;
  }
  
  // overwrite rho_i at the boundary
  const CFreal rho = this->_nodalPressure/RTovM;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    nstate[i] = _ys[i]*rho;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
