#include "NullDiffusiveEntity.hh"
#include "DiffusiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullDiffusiveEntity,
                       FiniteElementMethodData,
                       DiffusiveEntity,
                       FiniteElementModule>
NullDiffusiveEntityProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullDiffusiveEntity::NullDiffusiveEntity(const std::string& name) :
  DiffusiveEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullDiffusiveEntity::~NullDiffusiveEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& NullDiffusiveEntity::operator() ()
{
  _result = 0.;
  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

