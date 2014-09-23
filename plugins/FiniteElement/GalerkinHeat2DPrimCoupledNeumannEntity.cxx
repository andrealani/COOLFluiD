#include "GalerkinHeat2DPrimCoupledNeumannEntity.hh"
#include "NeumannEntity.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"
#include "Framework/SubSystemStatus.hh"
#include "Heat/Heat2DDiffusivePrim.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Heat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinHeat2DPrimCoupledNeumannEntity,
                       FiniteElementMethodData,
                       CoupledNeumannEntity,
                       FiniteElementModule>
GalerkinHeat2DDiffusivePrimCoupledNeumannEntity("GalerkinHeat2DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

GalerkinHeat2DPrimCoupledNeumannEntity::GalerkinHeat2DPrimCoupledNeumannEntity(const std::string& name) :
  CoupledNeumannEntity(name),
  _heatDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinHeat2DPrimCoupledNeumannEntity::~GalerkinHeat2DPrimCoupledNeumannEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinHeat2DPrimCoupledNeumannEntity::setup()
{
  CoupledNeumannEntity::setup();

  _heatDiffVarSet = _diffVarSet.d_castTo<Heat2DDiffusivePrim>();

}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinHeat2DPrimCoupledNeumannEntity::operator() ()
{
//   std::cout << FromHere().str() << std::endl;

  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
//   const std::vector<Framework::State*>& vars = *(_localElemData->solValues);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
//   const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

  RealVector& state = *((*_localElemData->solValues)[iQuadPoint]);
  RealVector& coord = *((*_localElemData->coord)[iQuadPoint]);

  const CFreal thickness =
    _heatDiffVarSet->getModel()->getThickness(*((*_localElemData->coord)[iQuadPoint]));
// CFout << "Thickness: " << thickness << "\n";
  cf_assert(_isAcceptedIdx < _isAccepted.size());
  if((_isAccepted[_isAcceptedIdx]) >= 0.)
  {
//     CFout << "Data: " << _index << "/" << _interfaceData.size() << "\n";
//     CFout << "Accepted: " << _isAcceptedIdx << "/" << _isAccepted.size()  << "\n";
//     CFout << "Coord: " << coord << "\n";
//     CFout << "State: " << state << "\n";
    if(_isRobin)
    {
      for(CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
        _varsRobin[iEq] = state[iEq];
      }

//       std::cout << "Data size: " << _interfaceData[_index].size() << "\n";
      for(CFuint iData = 0; iData < _interfaceData[_index].size(); ++iData) {
        _varsRobin[_nbEqs + iData ] = (_interfaceData[_index])[iData];
      }

      _vfRobin->evaluate(_varsRobin,_result);
//       std::cout << "Coord: " << coord << "\n";
//       std::cout << "Temperature: " << state << "\n";
//       std::cout << "Data: " << _interfaceData[_index] << "\n";
//       std::cout << "Flux: " << _result << "\n";
//       std::cout << std::endl;
      _index++;
    }
    else
    {
      _result = _interfaceData[_index];
      _index++;
    }
  }
  else
  {
    /// @note why not derive from NeumannEntity and if not accepted, use NeumannEntity
    //  CFout << "Refused: " << _isAcceptedIdx << "\n";
    if(_isRobin)
    {
      ///@todo this should do something by default
      throw Common::NotImplementedException(FromHere(),"Robin BC is not implemented for not accepted points");
    }
    else
    {
//       _result = 0.;
//       CFout << "Not accepted : " << state.getCoordinates() << "\n";
      // Variables are Coord, Time, States, Normals
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

