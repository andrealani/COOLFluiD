#include "GalerkinIndepSourceEntity.hh"
#include "IndepSourceEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealVector.hh"
#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinIndepSourceEntity,
                       FiniteElementMethodData,
                       IndepSourceEntity,
                       FiniteElementModule>
galerkinIndepSourceEntityProvider("Galerkin");

//////////////////////////////////////////////////////////////////////////////

GalerkinIndepSourceEntity::GalerkinIndepSourceEntity(const std::string& name) :
IndepSourceEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinIndepSourceEntity::~GalerkinIndepSourceEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinIndepSourceEntity::operator() ()
{
  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const Framework::State& vars = *((*_localElemData->solValues)[iQuadPoint]);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];

  _sourceVarSet->getIndepSourceCoefs(vars,_coefs);

  // this is assuming Galerkin approach,
  // and all shape functions are equal per equation
  //_result = shapeF[iState]*_coefs;
  _result = shapeF[iState]*_coefs;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////

