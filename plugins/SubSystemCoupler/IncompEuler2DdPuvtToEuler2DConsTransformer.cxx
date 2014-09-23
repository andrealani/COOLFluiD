#include "IncompEuler2DdPuvtToEuler2DConsTransformer.hh"
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


MethodStrategyProvider<IncompEuler2DdPuvtToEuler2DConsTransformer,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
incompEuler2DdPuvtToEuler2DConsTransformerProvider("IncompEuler2DdPuvtToEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtToEuler2DConsTransformer::IncompEuler2DdPuvtToEuler2DConsTransformer(const std::string& name) :
  PreVariableTransformer(name),
  _updateVar(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtToEuler2DConsTransformer::~IncompEuler2DdPuvtToEuler2DConsTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtToEuler2DConsTransformer::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  PreVariableTransformer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtToEuler2DConsTransformer::setup()
{
  CFAUTOTRACE;

  PreVariableTransformer::setup();

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  Common::SafePtr<SpaceMethodData> spaceMethodData = spaceMethod->getSpaceMethodData();
  _updateVar = (spaceMethodData->getUpdateVar()).d_castTo<Physics::NavierStokes::EulerVarSet>();

  //cf_assert(_updateVar->getName() == "dPuvt");

}
//////////////////////////////////////////////////////////////////////////////

RealVector* IncompEuler2DdPuvtToEuler2DConsTransformer::preTransform(const std::vector<GeoEntityIdx>& faces,const  RealVector& coord, const RealVector& original)
{
  cf_assert(original.size() == 4);
  _transVector.resize(4);

  const CFreal Cv = _updateVar->getModel()->getCv();
  const CFreal dP  = original[0];
  const CFreal U = original[1];
  const CFreal V = original[2];
  const CFreal T = original[3];
// unused // const CFreal pInf = _updateVar->getModel()->getThermodynamPressInf();
// unused //  const CFreal V2 = (U*U + V*V);
// unused //  const CFreal gamma = _updateVar->getModel()->getGamma();
// unused //  const CFreal R =_updateVar->getModel()->getR();
// unused //  const CFreal pressure = pInf + dP;

  const CFreal rho = _updateVar->getModel()->getDensity(dP,T);
  const CFreal rhoU = rho * U;
  const CFreal rhoV = rho * V;
  const CFreal rhoE = rho*(Cv*T + 0.5*(U*U + V*V));

  _transVector[0] = rho;
  _transVector[1] = rhoU;
  _transVector[2] = rhoV;
  _transVector[3] = rhoE;

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
