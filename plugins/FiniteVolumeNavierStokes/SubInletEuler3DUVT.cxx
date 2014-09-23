#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletEuler3DUVT.hh"
#include "Framework/MethodCommandProvider.hh"

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

MethodCommandProvider<SubInletEuler3DUVT, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletEuler3DUVTFVMCCProvider("SubInletEuler3DUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler3DUVT::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Vz","z velocity");
   options.addConfigOption< CFreal >("Vy","y velocity");
   options.addConfigOption< CFreal >("T","static temperature");
   options.addConfigOption< CFreal >("Vx","x velocity");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler3DUVT::SubInletEuler3DUVT(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _uinf = 0.0;
   setParameter("Vx",&_uinf);

  _vinf = 0.0;
   setParameter("Vy",&_vinf);

  _winf = 0.0;
   setParameter("Vz",&_winf);

  _temperature = 0.0;
   setParameter("T",&_temperature);
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler3DUVT::~SubInletEuler3DUVT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler3DUVT::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  // physical constants
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.0);
  const CFreal R = _varSet->getModel()->getR();
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];

  _dataGhostState[EulerTerm::RHO] = pInnerState/(R*_temperature);
  _dataGhostState[EulerTerm::VX] = 2.0*_uinf - _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = 2.0*_vinf - _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::VZ] = 2.0*_winf - _dataInnerState[EulerTerm::VZ];
  _dataGhostState[EulerTerm::P] = pInnerState;
  _dataGhostState[EulerTerm::V] = sqrt(_dataInnerState[EulerTerm::VX]*
				       _dataInnerState[EulerTerm::VX] +
				       _dataInnerState[EulerTerm::VY]*
				       _dataInnerState[EulerTerm::VY] +
				       _dataInnerState[EulerTerm::VZ]*
				       _dataInnerState[EulerTerm::VZ]);

  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
       _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
 
  _dataGhostState[EulerTerm::T] = 2.0*_temperature - _dataInnerState[EulerTerm::T];
  
  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler3DUVT::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();

  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _temperature/=_varSet->getModel()->getTempRef();
  _uinf /= _varSet->getModel()->getVelRef();
  _vinf /= _varSet->getModel()->getVelRef();
  _winf /= _varSet->getModel()->getVelRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
