#include "GalerkinConvectiveEntity.hh"
#include "ConvectiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinConvectiveEntity,
                       FiniteElementMethodData,
                       ConvectiveEntity,
                       FiniteElementModule>
galerkinConvectiveEntityProvider("Galerkin");

//////////////////////////////////////////////////////////////////////////////

GalerkinConvectiveEntity::GalerkinConvectiveEntity(const std::string& name) :
  ConvectiveEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinConvectiveEntity::~GalerkinConvectiveEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinConvectiveEntity::operator() ()
{
  /// @todo Very ineficient but we will refactor later when we improve the
  /// efficiency of the VolumeIntegrator

// unused //  const CFuint iState = _localElemData->iState;
// unused //  const CFuint jState = _localElemData->jState;

// unused //  const CFuint iQuadPoint =  _localElemData->quadPointID;
// unused //  const std::vector<Framework::State*>& vars = *(_localElemData->solValues);
// unused //  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
// unused //  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
// unused //  Framework::GeometricEntity* const geo = _localElemData->cell;

//  _convVarSet->getWeakConvMat(shapeF,shapeF,grad,grad,vars,geo,_result);

  throw Common::ShouldNotBeHereException (FromHere(),"Called GalerkinConvectiveEntity::operator()");

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

