#include "GalerkinStructMech3DDispInertiaEntity.hh"
#include "InertiaEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "FiniteElement/FiniteElementStructMech.hh"
#include "StructMech/StructMech3DInertiaDisp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMech;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMech3DDispInertiaEntity,
                       FiniteElementMethodData,
                       InertiaEntity,
                       FiniteElementStructMechModule>
GalerkinStructMech3DDispInertiaEntityProvider("GalerkinStructMech3DInertiaDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech3DDispInertiaEntity::GalerkinStructMech3DDispInertiaEntity(const std::string& name) :
  InertiaEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech3DDispInertiaEntity::~GalerkinStructMech3DDispInertiaEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech3DDispInertiaEntity::setup()
{
  InertiaEntity::setup();

  _structInertiaVarSet = _inertiaVarSet.d_castTo<StructMech3DInertiaDisp>();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinStructMech3DDispInertiaEntity::operator() ()
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
  _result(2,2) = shapeF[iState] * rho * shapeF[jState];
  _result(0,1) = 0;
  _result(0,2) = 0;
  _result(1,0) = 0;
  _result(1,2) = 0;
  _result(2,0) = 0;
  _result(2,1) = 0;

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////


