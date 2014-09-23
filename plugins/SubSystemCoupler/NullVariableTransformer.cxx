#include "Framework/PhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/NullVariableTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullVariableTransformer,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerModule>
nullVariableTransformerProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullVariableTransformer::NullVariableTransformer(const std::string& name) :
  PostVariableTransformer(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullVariableTransformer::~NullVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
