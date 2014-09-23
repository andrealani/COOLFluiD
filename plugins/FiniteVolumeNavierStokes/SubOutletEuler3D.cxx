#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEuler3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubOutletEuler3D, CellCenterFVMData, FiniteVolumeNavierStokesModule> subOutletEuler3DFVMCCProvider("SubOutletEuler3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler3D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler3D::SubOutletEuler3D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _pressure = 0.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler3D::~SubOutletEuler3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler3D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  _dataGhostState[EulerTerm::VX] = _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::VZ] = _dataInnerState[EulerTerm::VZ];
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::P] = 2.0*_pressure - _dataInnerState[EulerTerm::P];
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(_varSet->getModel()->getGamma()*
				       _dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = _dataInnerState[EulerTerm::T];

  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler3D::setup()
{
   FVMCC_BC::setup();
	
   _varSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
   _varSet->getModel()->resizePhysicalData(_dataInnerState);
   _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _pressure/=_varSet->getModel()->getPressRef();
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
