#include "GalerkinHeat3DPrimDiffEntity.hh"
#include "DiffusiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"
#include "Heat/Heat3DDiffusivePrim.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Heat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinHeat3DPrimDiffEntity,
                       FiniteElementMethodData,
                       DiffusiveEntity,
                       FiniteElementModule>
GalerkinHeat3DPrimDiffEntityProvider("GalerkinHeat3DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

GalerkinHeat3DPrimDiffEntity::GalerkinHeat3DPrimDiffEntity(const std::string& name) :
  DiffusiveEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinHeat3DPrimDiffEntity::~GalerkinHeat3DPrimDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinHeat3DPrimDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _heatDiffVarSet = _diffVarSet.d_castTo<Heat3DDiffusivePrim>();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinHeat3DPrimDiffEntity::operator() ()
{
  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];

  const CFreal k =
    _heatDiffVarSet->getModel()->getConductivity(
       *((*_localElemData->coord)[iQuadPoint]),
       *((*_localElemData->solValues)[iQuadPoint]));

  _result = grad(iState,XX)*grad(jState,XX)
          + grad(iState,YY)*grad(jState,YY)
          + grad(iState,ZZ)*grad(jState,ZZ);

  _result *= k;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

