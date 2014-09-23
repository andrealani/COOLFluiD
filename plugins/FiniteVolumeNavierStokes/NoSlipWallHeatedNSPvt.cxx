#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallHeatedNSPvt.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallHeatedNSPvt, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
NoSlipWallHeatedNS2DPuvtPuvtFVMCCProvider("NoSlipWallHeatedNS2DPuvtFVMCC");

MethodCommandProvider<NoSlipWallHeatedNSPvt, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
NoSlipWallHeatedNS3DPvtPuvtFVMCCProvider("NoSlipWallHeatedNS3DPvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallHeatedNSPvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("zWallVelocity","Z-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("HeatFlux","dT/dn");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallHeatedNSPvt::NoSlipWallHeatedNSPvt(const std::string& name) :
  FVMCC_BC(name),
  _diffusiveVarSet(CFNULL),
  _dummyGradients()
{
  addConfigOptionsTo(this);
  
  _xWallVelocity = 0.;
  setParameter("xWallVelocity",&_xWallVelocity);
  
  _yWallVelocity = 0.;
  setParameter("yWallVelocity",&_yWallVelocity);
  
  _zWallVelocity = 0.;
  setParameter("zWallVelocity",&_zWallVelocity);
  
  _heatFlux = 0.;
  setParameter("HeatFlux",&_heatFlux);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallHeatedNSPvt::~NoSlipWallHeatedNSPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallHeatedNSPvt::setup()
{
  FVMCC_BC::setup();
  
  Common::SafePtr<EulerVarSet> varSet =
    getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
  
  cf_assert(varSet.isNotNull());
  
  _xWallVelocity /= varSet->getModel()->getVelRef();
  _yWallVelocity /= varSet->getModel()->getVelRef();
  _zWallVelocity /= varSet->getModel()->getVelRef();
  
  //@todo adimensionalize the heat fluxes!!!!!!!
  
  _diffusiveVarSet = getMethodData().getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  cf_assert(_diffusiveVarSet.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallHeatedNSPvt::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // Set Tghost such that Q = lambda * (Tghost - Tin)/(distance_G_In)
  //-> set TG = Q*d/lambda + Tin
  
  const CFreal distance = MathFunctions::getDistance(innerState->getCoordinates(),ghostState->getCoordinates());
  _diffusiveVarSet->setWallDistance(distance/2.);
  const CFreal dynamicViscosity = _diffusiveVarSet->getDynViscosity(*innerState, _dummyGradients);
  
  const CFreal lambda = _diffusiveVarSet->getThermConductivity(*innerState,dynamicViscosity);
  
  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = 2.*_xWallVelocity - (*innerState)[1];
  (*ghostState)[2] = 2.*_yWallVelocity - (*innerState)[2];
  if (PhysicalModelStack::getActive()->getDim() == DIM_3D) {
    (*ghostState)[3] = 2.*_zWallVelocity - (*innerState)[3];
    (*ghostState)[4] = (_heatFlux * distance / lambda) + (*innerState)[4];
  }
  else {
    (*ghostState)[3] = (_heatFlux * distance / lambda) + (*innerState)[3];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
