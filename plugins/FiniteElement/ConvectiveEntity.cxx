#include "ConvectiveEntity.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

ConvectiveEntity::ConvectiveEntity(const std::string& name) :
  FEMIntegrableEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ConvectiveEntity::~ConvectiveEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

