#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSPvtSimple.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/MeshData.hh"
  
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSPvtSimple, CellCenterFVMData, FiniteVolumeNavierStokesModule> noSlipWallIsothermalNSPvtSimpleFVMCCProvider("NoSlipWallIsothermalNSPvtSimpleFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSPvtSimple::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("TWall","Wall temperature");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSPvtSimple::NoSlipWallIsothermalNSPvtSimple
(const std::string& name) :
  FVMCC_BC(name),
  _updateVarSet(CFNULL)
{
   addConfigOptionsTo(this);
  _wallTemp = 0.0;
   setParameter("TWall",&_wallTemp);
}
      
//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSPvtSimple::~NoSlipWallIsothermalNSPvtSimple()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSPvtSimple::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
  
  cf_assert(_wallTemp > 0.0);
  // adimensionalize the temperature
  _wallTemp /= _updateVarSet->getModel()->getTempRef();
}
      
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSPvtSimple::setGhostState(GeometricEntity *const face)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
      
  // here a fix is needed in order to have always ghostT > 0
  // if ghostT < 0  then the inner value is set
  State *const innerState = face->getState(0);
  
  const CFreal innerT = (dim == DIM_3D) ? (*innerState)[4] : (*innerState)[3];
  CFreal ghostT = 2.*_wallTemp - innerT;
  if (ghostT < 0.) {
    ghostT = _wallTemp;
  }
  cf_assert(ghostT > 0.);
  
  // reset the ghost node with the new position
  State *const ghostState = face->getState(1);
  
  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = -(*innerState)[1];
  (*ghostState)[2] = -(*innerState)[2];
  
  if (dim == DIM_3D) {
    (*ghostState)[3] = -(*innerState)[3];
    (*ghostState)[4] = ghostT;
  }
  else {
    (*ghostState)[3] = ghostT;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
