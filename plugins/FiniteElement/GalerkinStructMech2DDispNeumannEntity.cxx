#include "GalerkinStructMech2DDispNeumannEntity.hh"
#include "NeumannEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElementStructMech.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::StructMech;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

// MethodStrategyProvider<GalerkinStructMech2DDispNeumannEntity,
//                        FiniteElementMethodData,
//                        NeumannEntity,
//                        FiniteElementStructMechModule>
// GalerkinStructMech2DDiffusiveAxiDisp("GalerkinStructMech2DDiffusiveAxiDisp");

MethodStrategyProvider<GalerkinStructMech2DDispNeumannEntity,
                       FiniteElementMethodData,
                       NeumannEntity,
                       FiniteElementStructMechModule>
GalerkinStructMech2DDiffusiveDisp("GalerkinStructMech2DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DDispNeumannEntity::GalerkinStructMech2DDispNeumannEntity(const std::string& name) :
  NeumannEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMech2DDispNeumannEntity::~GalerkinStructMech2DDispNeumannEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMech2DDispNeumannEntity::setup()
{
  NeumannEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMech2DDiffusiveVarSet>();

}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinStructMech2DDispNeumannEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

///@todo change this for Axisymm
  const CFreal thickness = 1.;
//      _structDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

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

