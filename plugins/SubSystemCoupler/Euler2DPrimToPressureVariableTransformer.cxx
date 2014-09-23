#include "Euler2DPrimToPressureVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"
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

MethodStrategyProvider<Euler2DPrimToPressureVariableTransformer,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
euler2DPrimToPressureVariableTransformerProvider("Euler2DPrimToPressure");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToPressureVariableTransformer::Euler2DPrimToPressureVariableTransformer(const std::string& name) :
  PreVariableTransformer(name)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToPressureVariableTransformer::~Euler2DPrimToPressureVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrimToPressureVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  _transVector.resize(getTransformedSize(4));
}

//////////////////////////////////////////////////////////////////////////////

RealVector* Euler2DPrimToPressureVariableTransformer::preTransform(const std::vector<GeoEntityIdx>& faces,const  RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 4);

  _transVector[0] = original[3];

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
