#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMovePivt.hh"

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

MethodStrategyProvider<DistanceBasedExtrapolatorGMovePivt,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMovePivtProvider("DistanceBasedGMovePivt");

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMovePivt::DistanceBasedExtrapolatorGMovePivt
(const std::string& name) :
  DistanceBasedExtrapolatorGMoveRhoivt(name),
  _rhos(),
  _Rspecies(),
  _RTYovM()
{
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMovePivt::~DistanceBasedExtrapolatorGMovePivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMovePivt::setup()
{
  DistanceBasedExtrapolatorGMoveRhoivt::setup();
  
  const CFuint nbSpecies = _library->getNbSpecies();
  _rhos.resize(nbSpecies);
  
  _Rspecies.resize(nbSpecies);
  _library->setRiGas(_Rspecies);
  _RTYovM.resize(nbSpecies);
}
      
//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMovePivt::transform(const RealVector& in,
						   RealVector& out)
{
  // transform [p_i v T Tv_i] variables into [y_i v T Tv_i]
  const CFuint nbSpecies = _library->getNbSpecies();
  
  // pressure from p_i
  _pressure = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _pressure += in[ie];
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = nbSpecies + dim;
  const CFreal T = in[TID];
  out[TID] = T;
  
  CFuint tvID = TID + 1;
  const CFuint nbTv = _library->getNbTempVib();
  for  (CFuint i = 0; i < nbTv; ++i, ++tvID) {
    out[tvID] = in[tvID];
  }
  if (_library->getNbTe() == 1) {
    out[tvID] = in[tvID];
  }
  
  CFreal* tVec = (nbTv == 0) ? CFNULL : &(const_cast<RealVector&>(in).ptr()[TID+1]);
  const CFreal Te = _library->getTe(T,tVec);
  
  CFreal ovRho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Ti = (ie > 0) ? T : Te;
    _rhos[ie] = in[ie]/(_Rspecies[ie]*Ti);
    ovRho += _rhos[ie];
  }
  ovRho = 1./ovRho;
  
  // Set the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    out[ie] = _rhos[ie]*ovRho;
    _ys[ie] = out[ie];
  }
  
  // Set the speed components
  for  (CFuint i = 0; i < dim; ++i) {
    const CFuint iVar = nbSpecies + i;
    out[iVar] = in[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void DistanceBasedExtrapolatorGMovePivt::transformBack(RealVector& nstate)
{
  DataHandle<CFint> trsID = socket_trsID.getDataHandle();
  
  // transform [y_i v p] into [p_i v T Tv_i] variables, T and Tv_i
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
  CFreal* tVec = (_library->getNbTempVib() > 0) ? &nstate[TID+1] : CFNULL;
  const CFreal Te = _library->getTe(T,tVec);
  const CFreal rho = _library->density(Tdim, pdim, tVec)/refData[EulerTerm::RHO];
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Ti = (ie > 0) ? T : Te;
    nstate[ie] = max(0.0, rho*_ys[ie]*_Rspecies[ie]*Ti);
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

void DistanceBasedExtrapolatorGMovePivt::setSpeciesVariables(RealVector& nstate)
{
  // nstate is in [p_i v T Tv] variables
  // we want to recalculate p_i from _nodalPressure keeping mass fractions _ys constant 
  
  // we want to keep the same y_i and recalculate rho in function of p 
  const CFuint nbSpecies = _library->getNbSpecies();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = nbSpecies + dim;
  
  CFreal* tVec = (_library->getNbTempVib() > 0) ? &nstate[TID+1] : CFNULL;
  const CFreal T = nstate[TID];
  const CFreal Te = _library->getTe(T,tVec);
  
  CFreal ovRho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Ti = (ie > 0) ? T : Te;
    _rhos[ie] = nstate[ie]/(_Rspecies[ie]*Ti);
    ovRho += _rhos[ie];
  }
  ovRho = 1./ovRho;
  
  // Set the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = _rhos[ie]*ovRho;
    const CFreal Ti = (ie > 0) ? T : Te;
    _RTYovM[ie] = _Rspecies[ie]*_ys[ie]*Ti;
  }
  
  // overwrite rho_i at the boundary
  const CFreal rho = this->_nodalPressure/_RTYovM.sum();
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal Ti = (i > 0) ? T : Te;
    nstate[i] = rho*_ys[i]*_Rspecies[i]*Ti;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
