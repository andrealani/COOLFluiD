#include "DiffusiveEntity.hh"
#include "ComputeDiffusiveTerm.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/DiffusiveVarSet.hh"
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

MethodStrategyProvider<ComputeDiffusiveTerm, FiniteElementMethodData,
ComputeDiffusiveTerm, FiniteElementModule> diffusiveTermStrategyProvider("DiffusiveTerm");

//////////////////////////////////////////////////////////////////////////////

ComputeDiffusiveTerm::ComputeDiffusiveTerm(const std::string& name) :
  ComputeTerm<FiniteElementMethodData>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeDiffusiveTerm::~ComputeDiffusiveTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::configure ( Config::ConfigArgs& args )
{
  Framework::ComputeTerm<FiniteElementMethodData>::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::setup()
{
  ComputeTerm<FiniteElementMethodData>::setup();
  
  _diffusiveEntity = getMethodData().getDiffusiveEntity();
  cf_assert(_diffusiveEntity.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::computeTerm(GeometricEntity* const cell, RealMatrix& result)
{
  // compute Gradients
  getMethodData().getFEMVolumeIntegrator()->
    integrateFastGeneralFEMEntityOnGeoEnt<DiffusiveEntity>
      (*_diffusiveEntity, result);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD
