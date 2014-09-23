#include "GalerkinStructMechHeat3DIndepSourceEntity.hh"
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
using namespace COOLFluiD::Physics::StructMechHeat;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMechHeat3DIndepSourceEntity,
                       FiniteElementMethodData,
                       IndepSourceEntity,
                       FiniteElementModule>
GalerkinStructMechHeat2DDiffusiveAxiDispIndepSourceEntityProvider("GalerkinStructMechHeat3DDiffusiveAxiDisp");

MethodStrategyProvider<GalerkinStructMechHeat3DIndepSourceEntity,
                       FiniteElementMethodData,
                       IndepSourceEntity,
                       FiniteElementModule>
GalerkinStructMechHeat2DDiffusiveDispIndepSourceEntityProvider("GalerkinStructMechHeat3DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat3DIndepSourceEntity::GalerkinStructMechHeat3DIndepSourceEntity(const std::string& name) :
IndepSourceEntity(name),
_structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat3DIndepSourceEntity::~GalerkinStructMechHeat3DIndepSourceEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMechHeat3DIndepSourceEntity::setup()
{
  IndepSourceEntity::setup();

  _structDiffVarSet = getMethodData().getDiffusiveVar().d_castTo<StructMechHeat3DDiffusiveVarSet>();

}


//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinStructMechHeat3DIndepSourceEntity::operator() ()
{
  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const Framework::State& vars = *((*_localElemData->solValues)[iQuadPoint]);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];

  const CFreal dWdX = grad(iState,XX);
  const CFreal dWdY = grad(iState,YY);
  const CFreal dWdZ = grad(iState,ZZ);

  _sourceVarSet->getIndepSourceCoefs(vars,_coefs);

  // this is assuming Galerkin approach,
  // and all shape functions are equal per equation

  const CFreal youngModulus = _structDiffVarSet->getModel()->getYoung();
  const CFreal alpha = _structDiffVarSet->getModel()->getThermalExpansionCoef();
  const CFreal initialTemp = _structDiffVarSet->getModel()->getInitialTemp();
  const CFreal dT = vars[3] - initialTemp;
  const CFreal thickness =
    _structDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

  _result[0] = -alpha*youngModulus*dT*dWdX;
  _result[1] = -alpha*youngModulus*dT*dWdY;
  _result[2] = -alpha*youngModulus*dT*dWdZ;
  _result[3] = 0.;

  _result *= thickness;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////

