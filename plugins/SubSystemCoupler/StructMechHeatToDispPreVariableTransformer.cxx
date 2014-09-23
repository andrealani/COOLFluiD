#include "StructMechHeatToDispPreVariableTransformer.hh"
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

MethodStrategyProvider<StructMechHeatToDispPreVariableTransformer,
                       SubSysCouplerData,
                       PreVariableTransformer,
                       SubSystemCouplerHeatModule>
StructMechHeatToDispPreVariableTransformerProvider("StructMechHeatToDisp");

//////////////////////////////////////////////////////////////////////////////

StructMechHeatToDispPreVariableTransformer::StructMechHeatToDispPreVariableTransformer(const std::string& name) :
  PreVariableTransformer(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeatToDispPreVariableTransformer::~StructMechHeatToDispPreVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeatToDispPreVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  _transVector.resize(getTransformedSize(3));
}

//////////////////////////////////////////////////////////////////////////////

RealVector* StructMechHeatToDispPreVariableTransformer::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 3);


  _transVector[0] = original[0];
  _transVector[1] = original[1];

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
