#include "NullComputeIndepSourceTerm.hh"
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

MethodStrategyProvider<NullComputeIndepSourceTerm,
                       FiniteElementMethodData,
                       ComputeIndepSourceTerm,
                       FiniteElementModule>
nullIndepSourceTermStrategyProvider("NullIndepSourceTerm");

//////////////////////////////////////////////////////////////////////////////

NullComputeIndepSourceTerm::NullComputeIndepSourceTerm(const std::string& name) :
  ComputeIndepSourceTerm(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullComputeIndepSourceTerm::~NullComputeIndepSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

