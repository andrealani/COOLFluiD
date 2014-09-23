#include "GalerkinStructMech3DDispDiffEntity.hh"
#include "DiffusiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElementStructMech.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMech;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMech3DDispDiffEntity,
                       FiniteElementMethodData,
                       DiffusiveEntity,
                       FiniteElementStructMechModule>
GalerkinStructMech3DDispDiffEntityProvider("GalerkinStructMech3DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech3DDispDiffEntity::GalerkinStructMech3DDispDiffEntity(const std::string& name) :
  DiffusiveEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech3DDispDiffEntity::~GalerkinStructMech3DDispDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech3DDispDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMech3DDiffusiveDisp>();

  _isAnisotropic = _structDiffVarSet->getModel()->isAnisotropic();

  _meshMovement = _structDiffVarSet->isMeshMovement();
  _meshMovementMethod = _structDiffVarSet->getMeshMovementMethod();

  _isNonLinear = _structDiffVarSet->isNonLinear();

  _stiffness = _structDiffVarSet->getStiffnessMat();
  _lambda = _structDiffVarSet->getLambdaCoef();
  _mu = _structDiffVarSet->getMuCoef();

  resetStiffness();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinStructMech3DDispDiffEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
// unused //  const std::vector<Framework::State*>& vars = *(_localElemData->solValues);
// unused //  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

  if(_meshMovement){
    modifyStiffness(geo);
  }

  if(!_isAnisotropic){

    _result(0,0) = grad(iState,XX) * _lambda * grad(jState,XX) +
                  grad(iState,XX) * 2. * _mu * grad(jState,XX) +
                  grad(iState,YY) * _mu * grad(jState,YY) +
                  grad(iState,ZZ) * _mu * grad(jState,ZZ) ;

    _result(1,1) = grad(iState,YY) * _lambda * grad(jState,YY) +
                  grad(iState,YY) * 2. * _mu * grad(jState,YY) +
                  grad(iState,XX) * _mu * grad(jState,XX) +
                  grad(iState,ZZ) * _mu * grad(jState,ZZ) ;

    _result(2,2) = grad(iState,ZZ) * _lambda * grad(jState,ZZ) +
                  grad(iState,ZZ) * 2. * _mu * grad(jState,ZZ) +
                  grad(iState,XX) * _mu * grad(jState,XX) +
                  grad(iState,YY) * _mu * grad(jState,YY) ;

    _result(0,1) = grad(iState,XX) * _lambda * grad(jState,YY) +
                  grad(iState,YY) * _mu * grad(jState,XX) ;

    _result(0,2) = grad(iState,XX) * _lambda * grad(jState,ZZ) +
                  grad(iState,ZZ) * _mu * grad(jState,XX) ;

    _result(1,0) = grad(iState,YY) * _lambda * grad(jState,XX) +
                  grad(iState,XX) * _mu * grad(jState,YY) ;

    _result(1,2) = grad(iState,YY) * _lambda * grad(jState,ZZ) +
                  grad(iState,ZZ) * _mu * grad(jState,YY) ;

    _result(2,0) = grad(iState,ZZ) * _lambda * grad(jState,XX) +
                  grad(iState,XX) * _mu * grad(jState,ZZ) ;

    _result(2,1) = grad(iState,ZZ) * _lambda * grad(jState,YY) +
                  grad(iState,YY) * _mu * grad(jState,ZZ) ;


    if (_isNonLinear){
      /// Add NonLinear Contribution
      // for the moment, just isotropic...
      // Precompute some terms
      std::vector<State*>* nodalStates = geo->getStates();
      CFreal dNdxU(0);
      CFreal dNdyU(0);
      CFreal dNdzU(0);
      CFreal dNdxV(0);
      CFreal dNdyV(0);
      CFreal dNdzV(0);
      CFreal dNdxW(0);
      CFreal dNdyW(0);
      CFreal dNdzW(0);

      for (CFuint k=0;k<nodalStates->size();++k){
        dNdxU += grad(k,XX) * ((*((*nodalStates)[k]))[XX]);
        dNdyU += grad(k,YY) * ((*((*nodalStates)[k]))[XX]);
        dNdzU += grad(k,ZZ) * ((*((*nodalStates)[k]))[XX]);
        dNdxV += grad(k,XX) * ((*((*nodalStates)[k]))[YY]);
        dNdyV += grad(k,YY) * ((*((*nodalStates)[k]))[YY]);
        dNdzV += grad(k,ZZ) * ((*((*nodalStates)[k]))[YY]);
        dNdxW += grad(k,XX) * ((*((*nodalStates)[k]))[ZZ]);
        dNdyW += grad(k,YY) * ((*((*nodalStates)[k]))[ZZ]);
        dNdzW += grad(k,ZZ) * ((*((*nodalStates)[k]))[ZZ]);
      }

      _result(0,0) += 0.5 * grad(iState,XX) * _lambda * grad(jState,XX) * dNdxU +
                    0.5 * grad(iState,XX) * _lambda * grad(jState,YY) * dNdyU +
                    0.5 * grad(iState,XX) * _lambda * grad(jState,ZZ) * dNdzU +
                    grad(iState,XX) * _mu * grad(jState,XX) * dNdxU           +
                    grad(iState,YY) * _mu * grad(jState,XX) * dNdyU           +
                    grad(iState,ZZ) * _mu * grad(jState,XX) * dNdzU ;

      _result(0,1) += 0.5 * grad(iState,XX) * _lambda * grad(jState,XX) * dNdxV +
                    0.5 * grad(iState,XX) * _lambda * grad(jState,YY) * dNdyV +
                    0.5 * grad(iState,XX) * _lambda * grad(jState,ZZ) * dNdzV +
                    grad(iState,XX) * _mu * grad(jState,XX) * dNdxV           +
                    grad(iState,YY) * _mu * grad(jState,XX) * dNdyV           +
                    grad(iState,ZZ) * _mu * grad(jState,XX) * dNdzV ;

      _result(0,2) += 0.5 * grad(iState,XX) * _lambda * grad(jState,XX) * dNdxW +
                    0.5 * grad(iState,XX) * _lambda * grad(jState,YY) * dNdyW +
                    0.5 * grad(iState,XX) * _lambda * grad(jState,ZZ) * dNdzW +
                    grad(iState,XX) * _mu * grad(jState,XX) * dNdxW           +
                    grad(iState,YY) * _mu * grad(jState,XX) * dNdyW           +
                    grad(iState,ZZ) * _mu * grad(jState,XX) * dNdzW ;

      _result(1,0) += 0.5 * grad(iState,YY) * _lambda * grad(jState,XX) * dNdxU +
                    0.5 * grad(iState,YY) * _lambda * grad(jState,YY) * dNdyU +
                    0.5 * grad(iState,YY) * _lambda * grad(jState,ZZ) * dNdzU +
                    grad(iState,XX) * _mu * grad(jState,YY) * dNdxU           +
                    grad(iState,YY) * _mu * grad(jState,YY) * dNdyU           +
                    grad(iState,ZZ) * _mu * grad(jState,YY) * dNdzU ;

      _result(1,1) += 0.5 * grad(iState,YY) * _lambda * grad(jState,XX) * dNdxV +
                    0.5 * grad(iState,YY) * _lambda * grad(jState,YY) * dNdyV +
                    0.5 * grad(iState,YY) * _lambda * grad(jState,ZZ) * dNdzV +
                    grad(iState,XX) * _mu * grad(jState,YY) * dNdxV           +
                    grad(iState,YY) * _mu * grad(jState,YY) * dNdyV           +
                    grad(iState,ZZ) * _mu * grad(jState,YY) * dNdzV ;

      _result(1,2) += 0.5 * grad(iState,YY) * _lambda * grad(jState,XX) * dNdxW +
                    0.5 * grad(iState,YY) * _lambda * grad(jState,YY) * dNdyW +
                    0.5 * grad(iState,YY) * _lambda * grad(jState,ZZ) * dNdzW +
                    grad(iState,XX) * _mu * grad(jState,YY) * dNdxW           +
                    grad(iState,YY) * _mu * grad(jState,YY) * dNdyW           +
                    grad(iState,ZZ) * _mu * grad(jState,YY) * dNdzW ;

      _result(2,0) += 0.5 * grad(iState,ZZ) * _lambda * grad(jState,XX) * dNdxU +
                    0.5 * grad(iState,ZZ) * _lambda * grad(jState,YY) * dNdyU +
                    0.5 * grad(iState,ZZ) * _lambda * grad(jState,ZZ) * dNdzU +
                    grad(iState,XX) * _mu * grad(jState,ZZ) * dNdxU           +
                    grad(iState,YY) * _mu * grad(jState,ZZ) * dNdyU           +
                    grad(iState,ZZ) * _mu * grad(jState,ZZ) * dNdzU ;

      _result(2,1) += 0.5 * grad(iState,ZZ) * _lambda * grad(jState,XX) * dNdxV +
                    0.5 * grad(iState,ZZ) * _lambda * grad(jState,YY) * dNdyV +
                    0.5 * grad(iState,ZZ) * _lambda * grad(jState,ZZ) * dNdzV +
                    grad(iState,XX) * _mu * grad(jState,ZZ) * dNdxV           +
                    grad(iState,YY) * _mu * grad(jState,ZZ) * dNdyV           +
                    grad(iState,ZZ) * _mu * grad(jState,ZZ) * dNdzV ;

      _result(2,2) += 0.5 * grad(iState,ZZ) * _lambda * grad(jState,XX) * dNdxW +
                    0.5 * grad(iState,ZZ) * _lambda * grad(jState,YY) * dNdyW +
                    0.5 * grad(iState,ZZ) * _lambda * grad(jState,ZZ) * dNdzW +
                    grad(iState,XX) * _mu * grad(jState,ZZ) * dNdxW           +
                    grad(iState,YY) * _mu * grad(jState,ZZ) * dNdyW           +
                    grad(iState,ZZ) * _mu * grad(jState,ZZ) * dNdzW ;

    }
  }

  /// For Anisotropic Material
  /// Linear Elasticity (and Geometrical NonLinearities)
  else {

    const CFreal dWdX = grad(iState,XX);
    const CFreal dWdY = grad(iState,YY);
    const CFreal dWdZ = grad(iState,ZZ);
    const CFreal dNdX = grad(jState,XX);
    const CFreal dNdY = grad(jState,YY);
    const CFreal dNdZ = grad(jState,ZZ);

    _result(0,0) = (dWdX*c11 + dWdY*c41 + dWdZ*c61) * (dNdX) +
                  (dWdX*c14 + dWdY*c44 + dWdZ*c64) * (dNdY) +
                  (dWdX*c16 + dWdY*c46 + dWdZ*c66) * (dNdZ);

    _result(0,1) = (dWdX*c12 + dWdY*c42 + dWdZ*c62) * (dNdY) +
                  (dWdX*c14 + dWdY*c44 + dWdZ*c64) * (dNdX) +
                  (dWdX*c15 + dWdY*c45 + dWdZ*c65) * (dNdZ);

    _result(0,2) = (dWdX*c13 + dWdY*c43 + dWdZ*c63) * (dNdZ) +
                  (dWdX*c15 + dWdY*c45 + dWdZ*c65) * (dNdY) +
                  (dWdX*c16 + dWdY*c46 + dWdZ*c66) * (dNdX);

    _result(1,0) = (dWdX*c41 + dWdY*c21 + dWdZ*c51) * (dNdX) +
                  (dWdX*c44 + dWdY*c24 + dWdZ*c54) * (dNdY) +
                  (dWdX*c46 + dWdY*c26 + dWdZ*c56) * (dNdZ);

    _result(1,1) = (dWdX*c42 + dWdY*c22 + dWdZ*c52) * (dNdY) +
                  (dWdX*c44 + dWdY*c24 + dWdZ*c54) * (dNdX) +
                  (dWdX*c45 + dWdY*c25 + dWdZ*c55) * (dNdZ);

    _result(1,2) = (dWdX*c43 + dWdY*c23 + dWdZ*c53) * (dNdZ) +
                  (dWdX*c45 + dWdY*c25 + dWdZ*c55) * (dNdY) +
                  (dWdX*c46 + dWdY*c26 + dWdZ*c56) * (dNdX);

    _result(2,0) = (dWdX*c61 + dWdY*c51 + dWdZ*c31) * (dNdX) +
                  (dWdX*c64 + dWdY*c54 + dWdZ*c34) * (dNdY) +
                  (dWdX*c66 + dWdY*c56 + dWdZ*c36) * (dNdZ);

    _result(2,1) = (dWdX*c62 + dWdY*c52 + dWdZ*c32) * (dNdY) +
                  (dWdX*c64 + dWdY*c54 + dWdZ*c34) * (dNdX) +
                  (dWdX*c65 + dWdY*c55 + dWdZ*c35) * (dNdZ);

    _result(2,2) = (dWdX*c63 + dWdY*c53 + dWdZ*c33) * (dNdZ) +
                  (dWdX*c65 + dWdY*c55 + dWdZ*c35) * (dNdY) +
                  (dWdX*c66 + dWdY*c56 + dWdZ*c36) * (dNdX);

    if (_isNonLinear){
      /// Add NonLinear Contribution
      // Precompute some terms
      std::vector<State*>* nodalStates = geo->getStates();
      CFreal dNdxU(0);
      CFreal dNdyU(0);
      CFreal dNdzU(0);

      CFreal dNdxV(0);
      CFreal dNdyV(0);
      CFreal dNdzV(0);

      CFreal dNdxW(0);
      CFreal dNdyW(0);
      CFreal dNdzW(0);

      for (CFuint k=0;k<nodalStates->size();++k){
        dNdxU += grad(k,XX) * ((*((*nodalStates)[k]))[XX]);
        dNdyU += grad(k,YY) * ((*((*nodalStates)[k]))[XX]);
        dNdzU += grad(k,ZZ) * ((*((*nodalStates)[k]))[XX]);
        dNdxV += grad(k,XX) * ((*((*nodalStates)[k]))[YY]);
        dNdyV += grad(k,YY) * ((*((*nodalStates)[k]))[YY]);
        dNdzV += grad(k,ZZ) * ((*((*nodalStates)[k]))[YY]);
        dNdxW += grad(k,XX) * ((*((*nodalStates)[k]))[ZZ]);
        dNdyW += grad(k,YY) * ((*((*nodalStates)[k]))[ZZ]);
        dNdzW += grad(k,ZZ) * ((*((*nodalStates)[k]))[ZZ]);
      }

      _result(0,0) += (dWdX*c11 + dWdY*c41 + dWdZ*c61) * (0.5*dNdX*dNdxU) +
                    (dWdX*c12 + dWdY*c42 + dWdZ*c62) * (0.5*dNdY*dNdyU) +
                    (dWdX*c13 + dWdY*c43 + dWdZ*c63) * (0.5*dNdZ*dNdzU) +
                    (dWdX*c14 + dWdY*c44 + dWdZ*c64) * (0.5*(dNdX*dNdyU+dNdY*dNdxU)) +
                    (dWdX*c16 + dWdY*c46 + dWdZ*c66) * (0.5*(dNdX*dNdzU+dNdZ*dNdxU));

      _result(0,1) += (dWdX*c11 + dWdY*c41 + dWdZ*c61) * (0.5*dNdX*dNdxV) +
                    (dWdX*c12 + dWdY*c42 + dWdZ*c62) * (0.5*dNdY*dNdyV) +
                    (dWdX*c13 + dWdY*c43 + dWdZ*c63) * (0.5*dNdZ*dNdzV) +
                    (dWdX*c14 + dWdY*c44 + dWdZ*c64) * (0.5*(dNdX*dNdyV+dNdY*dNdxV)) +
                    (dWdX*c15 + dWdY*c45 + dWdZ*c65) * (0.5*(dNdY*dNdzV+dNdZ*dNdyV));

      _result(0,2) += (dWdX*c11 + dWdY*c41 + dWdZ*c61) * (0.5*dNdX*dNdxW) +
                    (dWdX*c12 + dWdY*c42 + dWdZ*c62) * (0.5*dNdY*dNdyW) +
                    (dWdX*c13 + dWdY*c43 + dWdZ*c63) * (0.5*dNdZ*dNdzW) +
                    (dWdX*c15 + dWdY*c45 + dWdZ*c65) * (0.5*(dNdY*dNdzW+dNdZ*dNdyW)) +
                    (dWdX*c16 + dWdY*c46 + dWdZ*c66) * (0.5*(dNdX*dNdzW+dNdZ*dNdxW));

      _result(1,0) += (dWdY*c21 + dWdX*c41 + dWdZ*c51) * (0.5*dNdX*dNdxU) +
                    (dWdY*c22 + dWdX*c42 + dWdZ*c52) * (0.5*dNdY*dNdyU) +
                    (dWdY*c23 + dWdX*c43 + dWdZ*c53) * (0.5*dNdZ*dNdzU) +
                    (dWdY*c24 + dWdX*c44 + dWdZ*c54) * (0.5*(dNdX*dNdyU+dNdY*dNdxU)) +
                    (dWdY*c26 + dWdX*c46 + dWdZ*c56) * (0.5*(dNdX*dNdzU+dNdZ*dNdxU));

      _result(1,1) += (dWdY*c21 + dWdX*c41 + dWdZ*c51) * (0.5*dNdX*dNdxV) +
                    (dWdY*c22 + dWdX*c42 + dWdZ*c52) * (0.5*dNdY*dNdyV) +
                    (dWdY*c23 + dWdX*c43 + dWdZ*c53) * (0.5*dNdZ*dNdzV) +
                    (dWdY*c24 + dWdX*c44 + dWdZ*c54) * (0.5*(dNdX*dNdyV+dNdY*dNdxV)) +
                    (dWdY*c25 + dWdX*c45 + dWdZ*c55) * (0.5*(dNdY*dNdzV+dNdZ*dNdyV));

      _result(1,2) += (dWdY*c21 + dWdX*c41 + dWdZ*c51) * (0.5*dNdX*dNdxW) +
                    (dWdY*c22 + dWdX*c42 + dWdZ*c52) * (0.5*dNdY*dNdyW) +
                    (dWdY*c23 + dWdX*c43 + dWdZ*c53) * (0.5*dNdZ*dNdzW) +
                    (dWdY*c25 + dWdX*c45 + dWdZ*c55) * (0.5*(dNdY*dNdzW+dNdZ*dNdyW)) +
                    (dWdY*c26 + dWdX*c46 + dWdZ*c56) * (0.5*(dNdX*dNdzW+dNdZ*dNdxW));

      _result(2,0) += (dWdZ*c31 + dWdY*c51 + dWdX*c61) * (0.5*dNdX*dNdxU) +
                    (dWdZ*c32 + dWdY*c52 + dWdX*c62) * (0.5*dNdY*dNdyU) +
                    (dWdZ*c33 + dWdY*c53 + dWdX*c63) * (0.5*dNdZ*dNdzU) +
                    (dWdZ*c34 + dWdY*c54 + dWdX*c64) * (0.5*(dNdX*dNdyU+dNdY*dNdxU)) +
                    (dWdZ*c36 + dWdY*c56 + dWdX*c66) * (0.5*(dNdX*dNdzU+dNdZ*dNdxU));

      _result(2,1) += (dWdZ*c31 + dWdY*c51 + dWdX*c61) * (0.5*dNdX*dNdxV) +
                    (dWdZ*c32 + dWdY*c52 + dWdX*c62) * (0.5*dNdY*dNdyV) +
                    (dWdZ*c33 + dWdY*c53 + dWdX*c63) * (0.5*dNdZ*dNdzV) +
                    (dWdZ*c34 + dWdY*c54 + dWdX*c64) * (0.5*(dNdX*dNdyV+dNdY*dNdxV)) +
                    (dWdZ*c35 + dWdY*c55 + dWdX*c65) * (0.5*(dNdY*dNdzV+dNdZ*dNdyV));

      _result(2,2) += (dWdZ*c31 + dWdY*c51 + dWdX*c61) * (0.5*dNdX*dNdxW) +
                    (dWdZ*c32 + dWdY*c52 + dWdX*c62) * (0.5*dNdY*dNdyW) +
                    (dWdZ*c33 + dWdY*c53 + dWdX*c63) * (0.5*dNdZ*dNdzW) +
                    (dWdZ*c35 + dWdY*c55 + dWdX*c65) * (0.5*(dNdY*dNdzW+dNdZ*dNdyW)) +
                    (dWdZ*c36 + dWdY*c56 + dWdX*c66) * (0.5*(dNdX*dNdzW+dNdZ*dNdxW));
    }

  }

  if(_meshMovement){
    resetStiffness();
  }

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech3DDispDiffEntity::modifyStiffness(Framework::GeometricEntity* const geo)
{
  cf_assert(false);

  CFreal factor = 1.;

  /// Method 1:
  /// Stiffness inversely proportional to the cell area/volume
  /// -> OK if mesh is not too uniform
  if(_meshMovementMethod == "VolumeBased"){
    CFreal volume = geo->computeVolume();
    cf_assert(volume > 0.);

    factor /= volume;
  }

  /// Method 2:
  /// Stiffness inversely proportional to the average distance to the wall
  if(_meshMovementMethod == "DistanceBased")
  {
    //Get the datahandle
    ///@todo change this, this is a dirty access to a datahandle...
    vector<Framework::State*>* geoStates = geo->getStates();
    std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
    std::string datahandleName = nsp + "_wallDistance";
    DataHandle<CFreal> wallDistance =
      MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

    //Compute Average Distance
    CFreal distance = 0.;
    for(CFuint iState = 0; iState < geoStates->size(); iState++)
    {
      distance += wallDistance[((*geoStates)[iState])->getLocalID()];
    }
    distance /= geoStates->size();
    cf_assert(distance > 0.);
    factor /= distance;
  }

  /// Method 3:
  /// Stiffness proportional to the cell quality
  if(_meshMovementMethod == "QualityBased")
  {
    //Get the datahandle
    ///@todo change this, this is a dirty access to a datahandle...
    std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
    std::string datahandleName = nsp + "_qualityCell";
    DataHandle<CFreal> cellQuality =
      MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

    //get the TRS index from the localID
    CFuint iGeo = geo->getID();
    CFreal quality = cellQuality[iGeo];

    cf_assert(quality > 1.);
    factor /= quality;
  }

  //Modify the stiffness
  _stiffness *= factor;

}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech3DDispDiffEntity::resetStiffness()
{

  /// elements of the constitutive matrix
  c11 = _stiffness(0,0);
  c12 = _stiffness(0,1);
  c13 = _stiffness(0,2);
  c14 = _stiffness(0,3);
  c15 = _stiffness(0,4);
  c16 = _stiffness(0,5);
  c21 = _stiffness(1,0);
  c22 = _stiffness(1,1);
  c23 = _stiffness(1,2);
  c24 = _stiffness(1,3);
  c25 = _stiffness(1,4);
  c26 = _stiffness(1,5);
  c31 = _stiffness(2,0);
  c32 = _stiffness(2,1);
  c33 = _stiffness(2,2);
  c34 = _stiffness(2,3);
  c35 = _stiffness(2,4);
  c36 = _stiffness(2,5);
  c41 = _stiffness(3,0);
  c42 = _stiffness(3,1);
  c43 = _stiffness(3,2);
  c44 = _stiffness(3,3);
  c45 = _stiffness(3,4);
  c46 = _stiffness(3,5);
  c51 = _stiffness(4,0);
  c52 = _stiffness(4,1);
  c53 = _stiffness(4,2);
  c54 = _stiffness(4,3);
  c55 = _stiffness(4,4);
  c56 = _stiffness(4,5);
  c61 = _stiffness(5,0);
  c62 = _stiffness(5,1);
  c63 = _stiffness(5,2);
  c64 = _stiffness(5,3);
  c65 = _stiffness(5,4);
  c66 = _stiffness(5,5);

  ///@todo reset the coeficients
  if(_meshMovement){
    cf_assert(false);
  }
//   _lambda =;
//   _mu =;

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

