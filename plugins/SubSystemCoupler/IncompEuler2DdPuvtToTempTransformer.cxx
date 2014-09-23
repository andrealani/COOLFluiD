#include "IncompEuler2DdPuvtToTempTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////


MethodStrategyProvider<IncompEuler2DdPuvtToTempTransformer,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
incompEuler2DdPuvtToTempTransformerProvider("IncompEuler2DdPuvtToTemp");

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtToTempTransformer::IncompEuler2DdPuvtToTempTransformer(const std::string& name) :
  PreVariableTransformer(name),
  _updateVar(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtToTempTransformer::~IncompEuler2DdPuvtToTempTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtToTempTransformer::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  PreVariableTransformer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtToTempTransformer::setup()
{
  CFAUTOTRACE;

  PreVariableTransformer::setup();

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  Common::SafePtr<SpaceMethodData> spaceMethodData = spaceMethod->getSpaceMethodData();
  _updateVar = (spaceMethodData->getUpdateVar()).d_castTo<Physics::NavierStokes::EulerVarSet>();

  //cf_assert(_updateVar->getName() == "dPuvt");
  _transVector.resize(1);
}
//////////////////////////////////////////////////////////////////////////////

RealVector* IncompEuler2DdPuvtToTempTransformer::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original)
{

  cf_assert(original.size() == 4);

  const CFreal T = original[3];

  _transVector[0] = T;

  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
