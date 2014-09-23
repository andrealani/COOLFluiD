#include "FiniteElement/DiffusiveEntity.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

DiffusiveEntity::DiffusiveEntity(const std::string& name) :
  FEMIntegrableEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffusiveEntity::~DiffusiveEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement
