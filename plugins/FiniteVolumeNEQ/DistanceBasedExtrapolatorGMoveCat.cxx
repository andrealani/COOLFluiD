#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveCat.hh"

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

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveCat,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveCatProvider("DistanceBasedGMoveCat");

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCat::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveCat::DistanceBasedExtrapolatorGMoveCat
(const std::string& name) :
  DistanceBasedExtrapolatorGMoveRhoivt(name)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveCat::~DistanceBasedExtrapolatorGMoveCat()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCat::applyBC()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle<CFint> trsID = socket_trsID.getDataHandle();
  
  const CFuint nbSpecies = _library->getNbSpecies();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = nbSpecies + dim;
  
  for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
    const CFuint idx = _wValuesIdx[i];
    if (idx < TID) {
      nodalStates[_currNodeID][idx] = _wValues[i];
    }
    else {
      // temperature is fixed only if radiative equilibrium is not in use
      if (!_radEquilibrium || _doNotApplyRadEq[_currNodeID]) {
	nodalStates[_currNodeID][idx] = (!_nodalValuesIDFlags[idx]) ?
	  _wValues[i] : getNodalValue(_orderedTrsList[trsID[_currNodeID]],idx,_currNodeID);
      }
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveCat::setup()
{
  DistanceBasedExtrapolatorGMoveRhoivt::setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
