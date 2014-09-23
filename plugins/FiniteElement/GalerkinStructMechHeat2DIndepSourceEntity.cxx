#include "GalerkinStructMechHeat2DIndepSourceEntity.hh"
#include "IndepSourceEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealVector.hh"
#include "FiniteElement/FiniteElementStructMechHeat.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMechHeat;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMechHeat2DIndepSourceEntity,
                       FiniteElementMethodData,
                       IndepSourceEntity,
                       FiniteElementStructMechHeatModule>
GalerkinHeatDeformationsIndepSourceEntityProvider("2DHeatDeformations");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DIndepSourceEntity::GalerkinStructMechHeat2DIndepSourceEntity(const std::string& name) :
IndepSourceEntity(name),
_structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DIndepSourceEntity::~GalerkinStructMechHeat2DIndepSourceEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMechHeat2DIndepSourceEntity::setup()
{
  CFAUTOTRACE;

  IndepSourceEntity::setup();

  _structDiffVarSet = getMethodData().getDiffusiveVar().d_castTo<StructMechHeat2DDiffusiveVarSet>();
  cf_assert(_structDiffVarSet.isNotNull());
}


//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinStructMechHeat2DIndepSourceEntity::operator() ()
{
  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const Framework::State& vars = *((*_localElemData->solValues)[iQuadPoint]);
//   const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];

  const CFreal dWdX = grad(iState,XX);
  const CFreal dWdY = grad(iState,YY);

  // this is assuming Galerkin approach,
  // and all shape functions are equal per equation
  const CFreal youngModulus = _structDiffVarSet->getModel()->getYoung();
  const CFreal poissonCoef = _structDiffVarSet->getModel()->getPoisson();
  const CFreal alpha = _structDiffVarSet->getModel()->getThermalExpansionCoef();
  const CFreal initialTemp = _structDiffVarSet->getModel()->getInitialTemp();
  const CFreal dT = vars[2] - initialTemp;
//   const CFreal radius = (*((*_localElemData->coord)[iQuadPoint]))[XX];
  const CFreal thickness =
    _structDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

  _result[0] = -alpha*youngModulus*dT*dWdX;
  _result[1] = -alpha*youngModulus*dT*dWdY;
  _result[2] = 0.;

  if(_structDiffVarSet->getModel()->isAxisymmetric())
  {
//if axi ???
//      const CFreal NoverR = shapeF[iState]/(radius+MathTools::MathConsts::CFrealEps());
//      _result[0] += alpha*youngModulus*dT*(NoverR);
    _result *= -1.;
  }

  _result *= thickness;
//if plane strain
  _result /= (1.+poissonCoef)*(1. - (2.*poissonCoef));

  if(_structDiffVarSet->getModel()->isAxisymmetric())
  {
    //if axi ???
    _result *= (1.+poissonCoef);
    //if plane strain
  }

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////

