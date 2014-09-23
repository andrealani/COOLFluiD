#include "GalerkinStructMechHeat2DDispCoupledNeumannEntity.hh"
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

MethodStrategyProvider<GalerkinStructMechHeat2DDispCoupledNeumannEntity,
                       FiniteElementMethodData,
                       CoupledNeumannEntity,
                       FiniteElementStructMechHeatModule>
GalerkinStructMechHeat2DDiffusiveAxiDispCoupled("GalerkinStructMechHeat2DDiffusiveAxiDisp");

MethodStrategyProvider<GalerkinStructMechHeat2DDispCoupledNeumannEntity,
                       FiniteElementMethodData,
                       CoupledNeumannEntity,
                       FiniteElementStructMechHeatModule>
GalerkinStructMechHeat2DDiffusiveDispCoupled("GalerkinStructMechHeat2DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DDispCoupledNeumannEntity::GalerkinStructMechHeat2DDispCoupledNeumannEntity(const std::string& name) :
  CoupledNeumannEntity(name),
  _structDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinStructMechHeat2DDispCoupledNeumannEntity::~GalerkinStructMechHeat2DDispCoupledNeumannEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinStructMechHeat2DDispCoupledNeumannEntity::setup()
{
  CoupledNeumannEntity::setup();

  _structDiffVarSet = _diffVarSet.d_castTo<StructMechHeat2DDiffusiveVarSet>();

}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinStructMechHeat2DDispCoupledNeumannEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
// unused //  const std::vector<Framework::State*>& vars = *(_localElemData->solValues);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
// unused //  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

  RealVector& state = *((*_localElemData->solValues)[iQuadPoint]);
  RealVector& coord = *((*_localElemData->coord)[iQuadPoint]);

  const CFreal thickness =
    _structDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));

  cf_assert(_isAcceptedIdx < _isAccepted.size());
  if((_isAccepted[_isAcceptedIdx]) >= 0.)
  {
/* CFout << "Data: " << _index << "/" << _interfaceData.size() << "\n";
 CFout << "Accepted: " << _isAcceptedIdx << "/" << _isAccepted.size()  << "\n";*/
 //CFout << "Coord: " << state.getCoordinates() << "\n";
    if(_isRobin)
    {

      for(CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
        _varsRobin[iEq] = state[iEq];
      }

      for(CFuint iData = 0; iData < _interfaceData[_index].size(); ++iData) {
        _varsRobin[_nbEqs + iData ] = (_interfaceData[_index])[iData];
      }

      _vfRobin->evaluate(_varsRobin,_result);
      _index++;
    }
    else{
      _result = _interfaceData[_index];
      _index++;
    }
  }
  else
  {
    /// why not derive from NeumannEntity and if not accepted, use NeumannEntity
    //  CFout << "Refused: " << _isAcceptedIdx << "\n";
    if(_isRobin)
    {
       _result = 0.;
    }
    else
    {
//       _result = 0.;
//       CFout << "Not accepted : " << state.getCoordinates() << "\n";
      ///Variables are Coord, Time, States, Normals
      // Set coordinates
      for(CFuint i = 0; i < _dim; ++i) {
        _vars[i] = coord[i];
      }

      // Set time
      _vars[_dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

      // Set state variables
      for(CFuint i = 0; i < _nbEqs; ++i) {
        _vars[i+_dim+1] = state[i];
      }

      // Set normals
      _normal = geo->computeAvgCellNormal();
      _normal.normalize();
      for(CFuint i = 0; i < _nbEqs; ++i) {
        _vars[i + _nbEqs + _dim + 1] = _normal[i];
      }

      _vf->evaluate(_vars,_result);
    }
  }
// CFout << "result * ShapeF[" << iState << "]\n";
  _result *= thickness * shapeF[iState];
  _isAcceptedIdx++;
  return _result;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

