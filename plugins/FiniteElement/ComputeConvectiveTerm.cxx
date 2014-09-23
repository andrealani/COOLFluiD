#include "ConvectiveEntity.hh"
#include "ComputeConvectiveTerm.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
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

MethodStrategyProvider<ComputeConvectiveTerm, FiniteElementMethodData, ComputeConvectiveTerm, FiniteElementModule> convectiveTermStrategyProvider("ConvectiveTerm");

//////////////////////////////////////////////////////////////////////////////

ComputeConvectiveTerm::ComputeConvectiveTerm(const std::string& name) :
  ComputeTerm<FiniteElementMethodData>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeConvectiveTerm::~ComputeConvectiveTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeConvectiveTerm::configure ( Config::ConfigArgs& args )
{
  Framework::ComputeTerm<FiniteElementMethodData>::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ComputeConvectiveTerm::computeTerm(GeometricEntity* const cell, RealMatrix& result)
{

  // compute Gradients
  getMethodData().getFEMVolumeIntegrator()->
    integrateFastGeneralFEMEntityOnGeoEnt<ConvectiveEntity>
      (*_convectiveEntity, result);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

