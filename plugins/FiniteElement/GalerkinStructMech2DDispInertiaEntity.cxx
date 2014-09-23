#include "GalerkinStructMech2DDispInertiaEntity.hh"
#include "InertiaEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "FiniteElement/FiniteElementStructMech.hh"
#include "StructMech/StructMech2DInertiaDisp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMech;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMech2DDispInertiaEntity,
                       FiniteElementMethodData,
                       InertiaEntity,
                       FiniteElementStructMechModule>
GalerkinStructMech2DDispInertiaEntityProvider("GalerkinStructMech2DInertiaDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DDispInertiaEntity::GalerkinStructMech2DDispInertiaEntity(const std::string& name) :
  InertiaEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DDispInertiaEntity::~GalerkinStructMech2DDispInertiaEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech2DDispInertiaEntity::setup()
{
  InertiaEntity::setup();

  _structInertiaVarSet = _inertiaVarSet.d_castTo<StructMech2DInertiaDisp>();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinStructMech2DDispInertiaEntity::operator() ()
{
  // this is assuming Galerkin approach,
  // and all shape functions are equal per equation
  const CFreal rho = _structInertiaVarSet->getModel()->getDensity();

  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
// unused //  const Framework::State& vars = *((*_localElemData->solValues)[iQuadPoint]);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];

  //2x2 matrix for the 2 equations of the problem
  //only diagonal entries are filled because no uv or vu terms
  //only u*u or v*v terms
  _result(0,0) = shapeF[iState] * rho * shapeF[jState];
  _result(1,1) = shapeF[iState] * rho * shapeF[jState];
  _result(0,1) = 0;
  _result(1,0) = 0;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////


