#include "ComputeLinearSourceTerm.hh"
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

MethodStrategyProvider<ComputeLinearSourceTerm, FiniteElementMethodData,
ComputeLinearSourceTerm, FiniteElementModule> linearSourceTermStrategyProvider("LinearSourceTerm");

//////////////////////////////////////////////////////////////////////////////

ComputeLinearSourceTerm::ComputeLinearSourceTerm(const std::string& name) :
  ComputeTerm<FiniteElementMethodData>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeLinearSourceTerm::~ComputeLinearSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeLinearSourceTerm::configure ( Config::ConfigArgs& args )
{
  Framework::ComputeTerm<FiniteElementMethodData>::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ComputeLinearSourceTerm::setup()
{
  CFAUTOTRACE;
  
  ComputeTerm<FiniteElementMethodData>::setup();
  
  _linearSourceEntity = getMethodData().getLinearSourceEntity();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeLinearSourceTerm::computeTerm(GeometricEntity* const cell, RealMatrix& result)
{
  getMethodData().getFEMVolumeIntegrator()->
    integrateFastGeneralFEMEntityOnGeoEnt<LinearSourceEntity>(*_linearSourceEntity, result);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
