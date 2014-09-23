#include "GalerkinStructMechHeat2DDispNeumannEntity.hh"
#include "NeumannEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElementStructMechHeat.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMechHeat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinStructMechHeat2DDispNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementStructMechHeatModule>
GalerkinStructMechHeat2DDiffusiveAxiDisp("GalerkinStructMechHeat2DDiffusiveAxiDisp");

MethodStrategyProvider<GalerkinStructMechHeat2DDispNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementStructMechHeatModule>
GalerkinStructMechHeat2DDiffusiveDisp("GalerkinStructMechHeat2DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DDispNeumannEntity::GalerkinStructMechHeat2DDispNeumannEntity(const std::string& name) :
  NeumannEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DDispNeumannEntity::~GalerkinStructMechHeat2DDispNeumannEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMechHeat2DDispNeumannEntity::setup()
{
  NeumannEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMechHeat2DDiffusiveVarSet>();

}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinStructMechHeat2DDispNeumannEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

  const CFreal thickness =
    _structDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

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

  _result *= thickness * shapeF[iState];
// CFout << "Coord: " << coord <<"\n";
// CFout << "Thickness: " << thickness <<"\n";
// CFout << "Result: " << _result <<"\n";

  return _result;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

