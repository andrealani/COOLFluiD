#include "NullComputeLinearSourceTerm.hh"
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

MethodStrategyProvider<NullComputeLinearSourceTerm,
                       FiniteElementMethodData,
                       ComputeLinearSourceTerm,
                       FiniteElementModule>
nullLinearSourceTermStrategyProvider("NullLinearSourceTerm");

//////////////////////////////////////////////////////////////////////////////

NullComputeLinearSourceTerm::NullComputeLinearSourceTerm(const std::string& name) :
  ComputeLinearSourceTerm(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullComputeLinearSourceTerm::~NullComputeLinearSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

