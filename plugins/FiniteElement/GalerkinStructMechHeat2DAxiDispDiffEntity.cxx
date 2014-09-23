#include "MathTools/MathConsts.hh"
#include "GalerkinStructMechHeat2DAxiDispDiffEntity.hh"
#include "DiffusiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElementStructMechHeat.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMechHeat;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMechHeat2DAxiDispDiffEntity,
                       FiniteElementMethodData,
                       DiffusiveEntity,
                       FiniteElementStructMechHeatModule>
GalerkinStructMechHeat2DAxiDispDiffEntityProvider("GalerkinStructMechHeat2DDiffusiveAxiDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DAxiDispDiffEntity::GalerkinStructMechHeat2DAxiDispDiffEntity(const std::string& name) :
  DiffusiveEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DAxiDispDiffEntity::~GalerkinStructMechHeat2DAxiDispDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMechHeat2DAxiDispDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMechHeat2DDiffusiveAxiDisp>();

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

RealMatrix& GalerkinStructMechHeat2DAxiDispDiffEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];

  const CFreal thickness =
    _structDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

  const CFreal radius = (*((*_localElemData->coord)[iQuadPoint]))[m_radiusID];

  const CFreal NoverR = shapeF[jState]/(radius + MathTools::MathConsts::CFrealEps());
  const CFreal WoverR = shapeF[iState]/(radius + MathTools::MathConsts::CFrealEps());
  const CFreal dWdR = grad(iState,m_radiusID);
  const CFreal dWdZ = grad(iState,m_zID);
  const CFreal dNdR = grad(jState,m_radiusID);
  const CFreal dNdZ = grad(jState,m_zID);

  /// For Isotropic Material
  /// Linear Elasticity
  if(!_isAnisotropic){
    _result(m_radiusID,m_radiusID) =    dWdR * _stiffness(0,0) * dNdR
                                      + dWdR * _stiffness(0,2) * NoverR
                                      + WoverR * _stiffness(2,0) * dNdR
                                      + WoverR * _stiffness(2,2) * NoverR
                                      + dWdZ * _stiffness(3,3) * dNdZ;

    _result(m_zID,m_radiusID) =    dWdZ * _stiffness(1,0) * dNdR
                                 + dWdZ * _stiffness(1,2) * NoverR
                                 + dWdR * _stiffness(3,3) * dNdZ;

    _result(m_radiusID,m_zID) =    dWdR * _stiffness(0,1) * dNdZ
                                 + WoverR * _stiffness(2,1) * dNdZ
                                 + dWdZ * _stiffness(3,3) * dNdR;

    _result(m_zID,m_zID) =    dWdZ * _stiffness(1,1) * dNdZ
                              + dWdR * _stiffness(3,3) * dNdR;

    if (_isNonLinear){
      cf_assert(false);
    }
  }

  /// For Anisotropic Material
  /// Linear Elasticity (and Geometrical NonLinearities)
  else{
    cf_assert(false);
  }

  const CFreal k = _structDiffVarSet->getModel()->getConductivity(
    *((*_localElemData->coord)[iQuadPoint]),*((*_localElemData->solValues)[iQuadPoint]));

  _result(2,2) = k * grad(iState,XX)*grad(jState,XX) + grad(iState,YY)*grad(jState,YY);

  _result(0,2) = 0.;
  _result(1,2) = 0.;
  _result(2,0) = 0.;
  _result(2,1) = 0.;

  _result *= thickness;
  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

