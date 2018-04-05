#include "ComputeIndepSourceTerm.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/SourceVarSet.hh"
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

MethodStrategyProvider<ComputeIndepSourceTerm, FiniteElementMethodData,
ComputeIndepSourceTerm, FiniteElementModule> indepSourceTermStrategyProvider("IndepSourceTerm");

//////////////////////////////////////////////////////////////////////////////

ComputeIndepSourceTerm::ComputeIndepSourceTerm(const std::string& name) :
  ComputeTerm<FiniteElementMethodData>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeIndepSourceTerm::~ComputeIndepSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeIndepSourceTerm::configure ( Config::ConfigArgs& args )
{
  Framework::ComputeTerm<FiniteElementMethodData>::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ComputeIndepSourceTerm::setup()
{
  ComputeTerm<FiniteElementMethodData>::setup();
  
 _indepSourceEntity = getMethodData().getIndepSourceEntity();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeIndepSourceTerm::computeTerm(GeometricEntity* const cell, RealVector& result)
{

  getMethodData().getFEMVolumeIntegrator()->
    integrateFastGeneralFEMEntityOnGeoEnt<IndepSourceEntity>(*_indepSourceEntity, result);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
