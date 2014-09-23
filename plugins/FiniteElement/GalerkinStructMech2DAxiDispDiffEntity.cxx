#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"
#include "FiniteElement/FiniteElementStructMech.hh"
#include "FiniteElement/GalerkinStructMech2DAxiDispDiffEntity.hh"
#include "FiniteElement/DiffusiveEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::StructMech;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMech2DAxiDispDiffEntity,
                       FiniteElementMethodData,
                       DiffusiveEntity,
                       FiniteElementStructMechModule>
GalerkinStructMech2DAxiDispDiffEntityProvider("GalerkinStructMech2DDiffusiveAxiDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DAxiDispDiffEntity::GalerkinStructMech2DAxiDispDiffEntity(const std::string& name) :
  DiffusiveEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DAxiDispDiffEntity::~GalerkinStructMech2DAxiDispDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech2DAxiDispDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMech2DDiffusiveAxiDisp>();

  _isAnisotropic = _structDiffVarSet->getModel()->isAnisotropic();

  _meshMovement = _structDiffVarSet->isMeshMovement();
  _meshMovementMethod = _structDiffVarSet->getMeshMovementMethod();

  _isNonLinear = _structDiffVarSet->isNonLinear();

  _stiffness = _structDiffVarSet->getStiffnessMat();

  std::string axisymmetryAxis = _structDiffVarSet->getAxisymmetryAxis();
  if((axisymmetryAxis == "X") || (axisymmetryAxis == "x")){
    m_radiusID = YY;
    m_zID = XX;
  }
  else{
    cf_assert((axisymmetryAxis == "Y") || (axisymmetryAxis == "y"));
    m_radiusID = XX;
    m_zID = YY;
  }

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinStructMech2DAxiDispDiffEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];

  const CFreal radius = (*((*_localElemData->coord)[iQuadPoint]))[m_radiusID];

  const CFreal NoverR = shapeF[iState]/(radius + MathTools::MathConsts::CFrealEps());
  const CFreal dWdR = grad(iState,m_radiusID);
  const CFreal dWdZ = grad(iState,m_zID);
  const CFreal dNdR = grad(jState,m_radiusID);
  const CFreal dNdZ = grad(jState,m_zID);

  /// For Isotropic Material
  /// Linear Elasticity
  if(!_isAnisotropic){
    _result(0,0) =    dWdR * _stiffness(0,0) * dNdR
                    + dWdR * _stiffness(0,2) * NoverR
                    + NoverR * _stiffness(2,0) * dNdR
                    + NoverR * _stiffness(2,2) * NoverR
                    + dWdZ * _stiffness(3,3) * dNdZ;

    _result(0,1) =    dWdZ * _stiffness(1,0) * dNdR
                    + dWdZ * _stiffness(1,2) * NoverR
                    + dWdR * _stiffness(3,3) * dNdZ;

    _result(1,0) =    dWdR * _stiffness(0,1) * dNdZ
                    + NoverR * _stiffness(2,1) * dNdZ
                    + dWdZ * _stiffness(3,3) * dNdR;

    _result(1,1) =    dWdR * _stiffness(0,0) * dNdR
                    + dWdR * _stiffness(0,2) * NoverR
                    + NoverR * _stiffness(2,0) * dNdR
                    + NoverR * _stiffness(2,2) * NoverR
                    + dWdZ * _stiffness(3,3) * dNdZ;

    if (_isNonLinear){
      cf_assert(false);
    }
  }

  /// For Anisotropic Material
  /// Linear Elasticity (and Geometrical NonLinearities)
  else{
    cf_assert(false);
  }

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

