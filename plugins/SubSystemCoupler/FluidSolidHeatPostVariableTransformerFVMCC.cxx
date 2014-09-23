#include "FluidSolidHeatPostVariableTransformerFVMCC.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethodData.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FluidSolidHeatPostVariableTransformerFVMCC,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerNavierStokesModule>
FluidSolidHeatPostVariableTransformerFVMCCProvider("FluidSolidHeatPreFVMCC");

//////////////////////////////////////////////////////////////////////

FluidSolidHeatPostVariableTransformerFVMCC::FluidSolidHeatPostVariableTransformerFVMCC(const std::string& name) :
  PostVariableTransformer(name),
  _varSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////

FluidSolidHeatPostVariableTransformerFVMCC::~FluidSolidHeatPostVariableTransformerFVMCC()
{
}

//////////////////////////////////////////////////////////////////////

void FluidSolidHeatPostVariableTransformerFVMCC::setup()
{

  PostVariableTransformer::setup();

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  cf_assert(spaceMethod.isNotNull());
  
  _varSet = spaceMethod->getSpaceMethodData()->getUpdateVar().d_castTo<EulerVarSet>();
  cf_assert(_varSet.isNotNull());
  
  _transVector.resize(1);
}


//////////////////////////////////////////////////////////////////////

RealVector* FluidSolidHeatPostVariableTransformerFVMCC::transform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original, const RealVector& pastTransformedVector)
{

  const CFreal tempRef = _varSet->getModel()->getTempRef();

  cf_assert(original.size() == 1);
  _transVector[0] = original[0]/tempRef;

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
