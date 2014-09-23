#include "StructMechHeatToTempPreVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "SubSystemCoupler/SubSystemCouplerHeat.hh"
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

MethodStrategyProvider<StructMechHeatToTempPreVariableTransformer,
                       SubSysCouplerData,
                       PreVariableTransformer,
                       SubSystemCouplerHeatModule>
StructMechHeatToTempPreVariableTransformerProvider("StructMechHeatToTemp");

//////////////////////////////////////////////////////////////////////////////

StructMechHeatToTempPreVariableTransformer::StructMechHeatToTempPreVariableTransformer(const std::string& name) :
  PreVariableTransformer(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeatToTempPreVariableTransformer::~StructMechHeatToTempPreVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeatToTempPreVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  _transVector.resize(getTransformedSize(3));
}

//////////////////////////////////////////////////////////////////////////////

RealVector* StructMechHeatToTempPreVariableTransformer::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 3);

  _transVector[0] = original[2];

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
