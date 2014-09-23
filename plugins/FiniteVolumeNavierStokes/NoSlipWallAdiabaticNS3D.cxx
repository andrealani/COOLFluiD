#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallAdiabaticNS3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler3DVarSet.hh"
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

MethodCommandProvider<NoSlipWallAdiabaticNS3D, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
noSlipWallAdiabaticNS3D3DFVMCCProvider("NoSlipWallAdiabaticNS3DFVMCC");
                                                                                                                                          
//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNS3D::NoSlipWallAdiabaticNS3D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNS3D::~NoSlipWallAdiabaticNS3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNS3D::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNS3D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSet->getModel()->getR();
  const CFreal ghostT = _dataInnerState[EulerTerm::P]/(R*_dataInnerState[EulerTerm::RHO]);
  const CFreal ghostP = _dataInnerState[EulerTerm::P];
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma / (gamma - 1.);

  _dataGhostState[EulerTerm::RHO] = ghostP/(R*ghostT);
  _dataGhostState[EulerTerm::VX] = -_dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = -_dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::VZ] = -_dataInnerState[EulerTerm::VZ];
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
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
