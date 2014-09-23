#include "GalerkinNeumannEntity.hh"
#include "NeumannEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementModule>
GalerkinStructMechHeat3DDispNeumannEntity("GalerkinStructMechHeat3DDiffusiveDisp");

MethodStrategyProvider<GalerkinNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementModule>
GalerkinStructMech3DDiffusiveDispNeumannEntity("GalerkinStructMech3DDiffusiveDisp");

MethodStrategyProvider<GalerkinNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementModule>
GalerkinHeat3DDiffusivePrimNeumannEntity("GalerkinHeat3DDiffusivePrim");

MethodStrategyProvider<GalerkinNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementModule>
GalerkinHeat2DDiffusivePrimNeumannEntity("GalerkinHeat2DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

GalerkinNeumannEntity::GalerkinNeumannEntity(const std::string& name) :
  NeumannEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinNeumannEntity::~GalerkinNeumannEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinNeumannEntity::setup()
{
  NeumannEntity::setup();

}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinNeumannEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
//   const std::vector<Framework::State*>& vars = *(_localElemData->solValues);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
//   const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

  ///Variables are Coord, Time, States, Normals
  // Set coordinates
  RealVector& coord = *((*_localElemData->coord)[iQuadPoint]);
  for(CFuint i = 0; i < _dim; ++i) {
    _vars[i] = coord[i];
  }

  // Set time
  _vars[_dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  // Set state variables
  RealVector& state = *((*_localElemData->solValues)[iQuadPoint]);
  for(CFuint i = 0; i < _nbEqs; ++i) {
    _vars[i+_dim+1] = state[i];
  }

  // Set normals
  ///@todo This is not correct, the direction of the face can be outward or inward!!!!
  _normal = geo->computeAvgCellNormal();
  _normal.normalize();
  for(CFuint i = 0; i < _dim; ++i) {
    _vars[i + _nbEqs + _dim + 1] = _normal[i];
  }

  _vf->evaluate(_vars,_result);

  _result *= shapeF[iState];

  return _result;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

