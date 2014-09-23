#include "GalerkinLinearSourceEntity.hh"
#include "LinearSourceEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathFunctions.hh"
#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinLinearSourceEntity,
                       FiniteElementMethodData,
                       LinearSourceEntity,
                       FiniteElementModule>
galerkinLinearSourceEntityProvider("Galerkin");

//////////////////////////////////////////////////////////////////////////////

GalerkinLinearSourceEntity::GalerkinLinearSourceEntity(const std::string& name) :
LinearSourceEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinLinearSourceEntity::~GalerkinLinearSourceEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinLinearSourceEntity::operator() ()
{
  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const Framework::State& vars = *((*_localElemData->solValues)[iQuadPoint]);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];

  _sourceVarSet->getLinearSourceCoefs(vars,_coefs);

  /// @todo is this correct?
  _result = shapeF[iState] * shapeF[jState] * _coefs;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////

