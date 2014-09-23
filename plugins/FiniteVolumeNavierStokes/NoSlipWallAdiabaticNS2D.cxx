#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallAdiabaticNS2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"

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

MethodCommandProvider< NoSlipWallAdiabaticNS2D,
                       CellCenterFVMData,
                       FiniteVolumeNavierStokesModule>
theNoSlipWallAdiabaticNS2D2DFVMCCProvider("NoSlipWallAdiabaticNS2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNS2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
   options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNS2D::NoSlipWallAdiabaticNS2D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  _xWallVelocity(0.),
  _yWallVelocity(0.)
{
   addConfigOptionsTo(this);
   setParameter("xWallVelocity",&_xWallVelocity);

   setParameter("yWallVelocity",&_yWallVelocity);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNS2D::~NoSlipWallAdiabaticNS2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNS2D::setup()
{

  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _xWallVelocity /= _varSet->getModel()->getVelRef();
  _yWallVelocity /= _varSet->getModel()->getVelRef();

}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNS2D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSet->getModel()->getR();
  const CFreal ghostT = _dataInnerState[EulerTerm::P]/
    (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal ghostP = _dataInnerState[EulerTerm::P];
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.);

  _dataGhostState[EulerTerm::RHO] = ghostP/(R*ghostT);
  _dataGhostState[EulerTerm::VX] = 2*_xWallVelocity - _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = 2*_yWallVelocity - _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				       _dataGhostState[EulerTerm::VX] +
				       _dataGhostState[EulerTerm::VY]*
				       _dataGhostState[EulerTerm::VY]);
  _dataGhostState[EulerTerm::P] = ghostP;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*ghostP +
			 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
    _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*ghostP/_dataGhostState[EulerTerm::RHO]);
   
  _dataGhostState[EulerTerm::T] = _dataInnerState[EulerTerm::T];

  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
