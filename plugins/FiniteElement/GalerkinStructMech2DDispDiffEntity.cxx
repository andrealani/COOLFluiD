#include "GalerkinStructMech2DDispDiffEntity.hh"
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

MethodStrategyProvider<GalerkinStructMech2DDispDiffEntity,
                       FiniteElementMethodData,
                            DiffusiveEntity,
                            FiniteElementStructMechModule>
GalerkinStructMech2DDispDiffEntityProvider("GalerkinStructMech2DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DDispDiffEntity::GalerkinStructMech2DDispDiffEntity(const std::string& name) :
  DiffusiveEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DDispDiffEntity::~GalerkinStructMech2DDispDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech2DDispDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMech2DDiffusiveDisp>();

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

RealMatrix& GalerkinStructMech2DDispDiffEntity::operator() ()
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

  const CFreal dWdX = grad(iState,XX);
  const CFreal dWdY = grad(iState,YY);
  const CFreal dNdX = grad(jState,XX);
  const CFreal dNdY = grad(jState,YY);

  /// For Isotropic Material
  /// Linear Elasticity (and Geometrical NonLinearities)
  if(!_isAnisotropic){
    _result(0,0) = dWdX * _c11 * dNdX + dWdY * _c66 * dNdY;
    _result(0,1) = dWdX * _c12 * dNdY + dWdY * _c66 * dNdX;
    _result(1,0) = dWdY * _c21 * dNdX + dWdX * _c66 * dNdY;
    _result(1,1) = dWdX * _c66 * dNdX + dWdY * _c22 * dNdY;

    if (_isNonLinear){
      /// Add NonLinear Contribution
      // Precompute some terms
      std::vector<State*>* nodalStates = geo->getStates();
      CFreal dNdxU(0);
      CFreal dNdyU(0);
      CFreal dNdxV(0);
      CFreal dNdyV(0);

      for (CFuint k=0;k<nodalStates->size();++k){
        dNdxU += grad(k,XX) * ((*((*nodalStates)[k]))[XX]);
        dNdxV += grad(k,XX) * ((*((*nodalStates)[k]))[YY]);
        dNdyU += grad(k,YY) * ((*((*nodalStates)[k]))[XX]);
        dNdyV += grad(k,YY) * ((*((*nodalStates)[k]))[YY]);
      }

      _result(0,0) += 0.5 * dWdX * _lambda * dNdX * dNdxU +
                     0.5 * dWdX * _lambda * dNdY * dNdyU +
                     dWdX * _mu * dNdX * dNdxU           +
                     dWdY * _mu * dNdX * dNdyU ;

      _result(0,1) += 0.5 * dWdX * _lambda * dNdX * dNdxV +
                     0.5 * dWdX * _lambda * dNdY * dNdyV +
                     dWdX * _mu * dNdX * dNdxV           +
                     dWdY * _mu * dNdX * dNdyV ;

      _result(1,0) += 0.5 * dWdY * _lambda * dNdY * dNdyU +
                     0.5 * dWdY * _lambda * dNdX * dNdxU +
                     dWdY * _mu * dNdY * dNdyU           +
                     dWdX * _mu * dNdY * dNdxU ;

      _result(1,1) += 0.5 * dWdY * _lambda * dNdY * dNdyV +
                     0.5 * dWdY * _lambda * dNdX * dNdxV +
                     dWdY * _mu * dNdY * dNdyV           +
                     dWdX * _mu * dNdY * dNdxV ;
    }
  }

  /// For Anisotropic Material
  /// Linear Elasticity (and Geometrical NonLinearities)
  else{
    _result(0,0) = (dWdX*_c11 + dWdY*_c61) * (dNdX) +
                  (dWdX*_c16 + dWdY*_c66) * (dNdY);

    _result(0,1) = (dWdX*_c12 + dWdY*_c62) * (dNdY) +
                  (dWdX*_c16 + dWdY*_c66) * (dNdX);

    _result(1,0) = (dWdX*_c61 + dWdY*_c21) * (dNdX) +
                  (dWdX*_c66 + dWdY*_c26) * (dNdY);

    _result(1,1) = (dWdX*_c62 + dWdY*_c22) * (dNdY) +
                  (dWdX*_c66 + dWdY*_c26) * (dNdX);

    if (_isNonLinear){
      /// Add NonLinear Contribution
      // Precompute some terms
      std::vector<State*>* nodalStates = geo->getStates();
      CFreal dNdxU(0);
      CFreal dNdyU(0);

      CFreal dNdxV(0);
      CFreal dNdyV(0);

      for (CFuint k=0;k<nodalStates->size();++k){
        dNdxU += grad(k,XX) * ((*((*nodalStates)[k]))[XX]);
        dNdyU += grad(k,YY) * ((*((*nodalStates)[k]))[XX]);
        dNdxV += grad(k,XX) * ((*((*nodalStates)[k]))[YY]);
        dNdyV += grad(k,YY) * ((*((*nodalStates)[k]))[YY]);
      }

      _result(0,0) += (dWdX*_c11 + dWdY*_c61) * (0.5*dNdX*dNdxU) +
                    (dWdX*_c12 + dWdY*_c62) * (0.5*dNdY*dNdyU) +
                    (dWdX*_c16 + dWdY*_c66) * (0.5*dNdX*dNdyU + 0.5*dNdY*dNdxU);

      _result(0,1) += (dWdX*_c11 + dWdY*_c61) * (0.5*dNdX*dNdxV) +
                    (dWdX*_c12 + dWdY*_c62) * (0.5*dNdY*dNdyV) +
                    (dWdX*_c16 + dWdY*_c66) * (0.5*dNdX*dNdyV + 0.5*dNdY*dNdxV);

      _result(1,0) += (dWdY*_c21 + dWdX*_c61) * (0.5*dNdX*dNdxU) +
                    (dWdY*_c22 + dWdX*_c62) * (0.5*dNdY*dNdyU) +
                    (dWdY*_c26 + dWdX*_c66) * (0.5*dNdX*dNdyU + 0.5*dNdY*dNdxU);

      _result(1,1) += (dWdY*_c21 + dWdX*_c61) * (0.5*dNdX*dNdxV) +
                    (dWdY*_c22 + dWdX*_c62) * (0.5*dNdY*dNdyV) +
                    (dWdY*_c26 + dWdX*_c66) * (0.5*dNdX*dNdyV + 0.5*dNdY*dNdxV);
    }
  }

  if(_meshMovement){
    resetStiffness();
  }

  return _result;
}


//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech2DDispDiffEntity::modifyStiffness(Framework::GeometricEntity* const geo)
{

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
  // only _c11,_c12... are used in the getWeakMat() so, dont change the stiffness RealMatrix but just these CFreal's
  _c11 *= factor;
  _c12 *= factor;
  _c16 *= factor;
  _c21 *= factor;
  _c22 *= factor;
  _c26 *= factor;
  _c61 *= factor;
  _c62 *= factor;
  _c66 *= factor;

}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech2DDispDiffEntity::resetStiffness()
{

  _c11 = _stiffness(0,0);
  _c12 = _stiffness(0,1);
  _c16 = _stiffness(0,2);
  _c21 = _stiffness(1,0);
  _c22 = _stiffness(1,1);
  _c26 = _stiffness(1,2);
  _c61 = _stiffness(2,0);
  _c62 = _stiffness(2,1);
  _c66 = _stiffness(2,2);

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

