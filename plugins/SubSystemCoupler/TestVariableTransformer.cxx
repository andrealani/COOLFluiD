#include "TestVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
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

MethodStrategyProvider<TestVariableTransformer,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerModule>
testVariableTransformerProvider("Test");

//////////////////////////////////////////////////////////////////////////////

TestVariableTransformer::TestVariableTransformer(const std::string& name) :
  PostVariableTransformer(name)
{
}

//////////////////////////////////////////////////////////////////////////////

TestVariableTransformer::~TestVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector* TestVariableTransformer::transform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original, const RealVector& pastTransformedVector)
{
  _transVector.resize(original.size());
  for(CFuint i=0;i<original.size();++i)
  {
    _transVector[i] = 10.*original[i];
  }
  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
