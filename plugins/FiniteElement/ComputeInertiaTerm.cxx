#include "InertiaEntity.hh"
#include "ComputeInertiaTerm.hh"
#include "Framework/VolumeIntegrator.hh"
#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ComputeInertiaTerm, FiniteElementMethodData,
ComputeInertiaTerm, FiniteElementModule> inertiaTermStrategyProvider("InertiaTerm");

//////////////////////////////////////////////////////////////////////////////

ComputeInertiaTerm::ComputeInertiaTerm(const std::string& name) :
  ComputeTerm<FiniteElementMethodData>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeInertiaTerm::~ComputeInertiaTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInertiaTerm::configure ( Config::ConfigArgs& args )
{
  Framework::ComputeTerm<FiniteElementMethodData>::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ComputeInertiaTerm::setup()
{

  _inertiaEntity = getMethodData().getInertiaEntity();
  cf_assert(_inertiaEntity.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void ComputeInertiaTerm::computeTerm(GeometricEntity* const cell, RealMatrix& result)
{

  getMethodData().getFEMVolumeIntegrator()->
    integrateFastGeneralFEMEntityOnGeoEnt<InertiaEntity>
      (*_inertiaEntity, result);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
