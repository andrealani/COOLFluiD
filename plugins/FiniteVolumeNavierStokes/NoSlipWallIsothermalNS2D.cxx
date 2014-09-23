#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNS2D.hh"
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

MethodCommandProvider<NoSlipWallIsothermalNS2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> noSlipWallIsothermalNS2DFVMCCProvider("NoSlipWallIsothermalNS2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNS2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("TWall","Wall temperature");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNS2D::NoSlipWallIsothermalNS2D(const std::string& name) :
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

NoSlipWallIsothermalNS2D::~NoSlipWallIsothermalNS2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNS2D::setup()
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

void NoSlipWallIsothermalNS2D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const Node& innerNode = innerState->getCoordinates();
  // store the ghost node
  _tempGhostNode = ghostState->getCoordinates();

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSet->getModel()->getR();
  const CFreal innerT = _dataInnerState[EulerTerm::P]/
    (R*_dataInnerState[EulerTerm::RHO]);

  // this middle node is by construction on the boundary face
  _midNode = 0.5*(innerNode + _tempGhostNode);

  // here a fix is needed in order to have always ghostT > 0
  // dynamic relocation of the ghost state: the position of the
  // ghost state is locally changed, and the BC is imposed
  // using a weighted average of ghost state (in thenew location)
  // and inner state
  CFreal ghostT = 2.*_wallTemp - innerT;
  CFreal factor = 1.;
  while (ghostT < 0.) {
    _tempNode = 0.5*(_midNode + _tempGhostNode);
    factor *= 2.;
    ghostT = ((factor + 1.)*_wallTemp - innerT)/factor;
    // move the ghost to the new position
    _tempGhostNode = _tempNode;
  }
  cf_assert(ghostT > 0.);

  // reset the ghost state
  ghostState->getCoordinates() = _tempGhostNode;

  const CFreal ghostP = _dataInnerState[EulerTerm::P];
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.);

  _dataGhostState[EulerTerm::RHO] = ghostP/(R*ghostT);
  _dataGhostState[EulerTerm::VX] = -_dataInnerState[EulerTerm::VX]/factor;
  _dataGhostState[EulerTerm::VY] = -_dataInnerState[EulerTerm::VY]/factor;
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V]/factor;
  _dataGhostState[EulerTerm::P] = ghostP;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*ghostP +
				   0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
    _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*ghostP/_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::T] = ghostT;
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
