#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermNS2D.hh"
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

MethodCommandProvider<NoSlipWallIsothermNS2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
noSlipWallIsothermNS2DFVMCCProvider("NoSlipWallIsothermNS2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermNS2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("TWall","Wall temperature");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermNS2D::NoSlipWallIsothermNS2D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  _tempNode(),
  _midNode(),
  _tempGhostNode()
{
   addConfigOptionsTo(this);
  _wallTemp = 0.0;
   setParameter("TWall",&_wallTemp);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermNS2D::~NoSlipWallIsothermNS2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermNS2D::setup()
{

  FVMCC_BC::setup();

  _tempNode.resize(PhysicalModelStack::getActive()->getDim());
  _midNode.resize(PhysicalModelStack::getActive()->getDim());
  _tempGhostNode.resize(PhysicalModelStack::getActive()->getDim());

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  cf_assert(_wallTemp > 0.0);
  // adimensionalize the temperature
  _wallTemp /= _varSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermNS2D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSet->getModel()->getR();
  const CFreal innerT = _dataInnerState[EulerTerm::P]/
    (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal ghostT = 2.*_wallTemp - innerT;
  const CFreal ghostP = _dataInnerState[EulerTerm::P];
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 =  gamma/(gamma - 1.);

  _dataGhostState[EulerTerm::RHO] = ghostP/(R*ghostT);
  _dataGhostState[EulerTerm::VX] = -_dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = -_dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::P] = ghostP;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*ghostP +
			 0.5*_dataGhostState[EulerTerm::RHO]*
			 _dataGhostState[EulerTerm::V]*
			 _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*ghostP/_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::T] = ghostT;
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermNS2D::computeFlux(RealVector& result)
{
  getMethodData().setUseAverageFlux(true);
  getMethodData().getFluxSplitter()->computeFlux(result);
  getMethodData().setUseAverageFlux(false);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
