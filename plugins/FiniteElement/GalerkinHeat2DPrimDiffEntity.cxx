#include "GalerkinHeat2DPrimDiffEntity.hh"
#include "DiffusiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"
#include "Heat/Heat2DDiffusivePrim.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Heat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinHeat2DPrimDiffEntity,
                       FiniteElementMethodData,
                       DiffusiveEntity,
                       FiniteElementModule>
GalerkinHeat2DPrimDiffEntityProvider("GalerkinHeat2DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

GalerkinHeat2DPrimDiffEntity::GalerkinHeat2DPrimDiffEntity(const std::string& name) :
  DiffusiveEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinHeat2DPrimDiffEntity::~GalerkinHeat2DPrimDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinHeat2DPrimDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _heatDiffVarSet = _diffVarSet.d_castTo<Heat2DDiffusivePrim>();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinHeat2DPrimDiffEntity::operator() ()
{

  cf_assert(_localElemData.isNotNull());
  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
  const CFreal k = _heatDiffVarSet->getModel()->getConductivity(
    *((*_localElemData->coord)[iQuadPoint]),*((*_localElemData->solValues)[iQuadPoint]));

  const CFreal thickness =
    _heatDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

  _result = grad(iState,XX)*grad(jState,XX) + grad(iState,YY)*grad(jState,YY);

  _result *= k * thickness;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

