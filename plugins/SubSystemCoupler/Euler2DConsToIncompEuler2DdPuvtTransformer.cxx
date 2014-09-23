#include "Euler2DConsToIncompEuler2DdPuvtTransformer.hh"
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


MethodStrategyProvider<Euler2DConsToIncompEuler2DdPuvtTransformer,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerNavierStokesModule>
Euler2DConsToIncompEuler2DdPuvtTransformerProvider("Euler2DConsToIncompEuler2DdPuvt");

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToIncompEuler2DdPuvtTransformer::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("p0Inf","Thermodynamic pressure infinity value.");
   options.addConfigOption< CFreal >("RDim","Constant of perfect gas.");
   options.addConfigOption< CFreal >("gamma","Specific heat ratio.");
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToIncompEuler2DdPuvtTransformer::Euler2DConsToIncompEuler2DdPuvtTransformer(const std::string& name) :
  PostVariableTransformer(name),
  _updateVar(CFNULL)
{
   addConfigOptionsTo(this);

   _p0Inf = 0.;
   setParameter("p0Inf",&_p0Inf);

   _RDim = 287.046;
   setParameter("RDim",&_RDim);

  _gamma = 1.4;
   setParameter("gamma",&_gamma);

}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToIncompEuler2DdPuvtTransformer::~Euler2DConsToIncompEuler2DdPuvtTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToIncompEuler2DdPuvtTransformer::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  PostVariableTransformer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToIncompEuler2DdPuvtTransformer::setup()
{
  CFAUTOTRACE;

  PostVariableTransformer::setup();

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  Common::SafePtr<SpaceMethodData> spaceMethodData = spaceMethod->getSpaceMethodData();
  _updateVar = (spaceMethodData->getUpdateVar()).d_castTo<Physics::NavierStokes::EulerVarSet>();

  _transVector.resize(getTransformedSize(4));

}
//////////////////////////////////////////////////////////////////////////////

RealVector* Euler2DConsToIncompEuler2DdPuvtTransformer::transform(const std::vector<GeoEntityIdx>& faces,const  RealVector& coord, const RealVector& original, const RealVector& pastTransformedVector)
{
  cf_assert(original.size() == 4);

  const CFreal rho = original[0];
  const CFreal rhoU = original[1];
  const CFreal rhoV = original[2];
  const CFreal rhoE = original[3];

  const CFreal rhoV2 = (rhoU*rhoU + rhoV*rhoV)/rho;
  const CFreal gamma = _updateVar->getModel()->getGamma();
  const CFreal R =_updateVar->getModel()->getR();

  const CFreal pressure = (gamma - 1.)*(rhoE - (0.5*rhoV2));
  const CFreal dP = pressure - _p0Inf;
  const CFreal temperature = pressure/(R*rho);

  _transVector[0] = dP;
  _transVector[1] = rhoU/rho;
  _transVector[2] = rhoV/rho;
  _transVector[3] = temperature;

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
