#include "GalerkinCoupledNeumannEntity.hh"
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

MethodStrategyProvider<GalerkinCoupledNeumannEntity,
                       FiniteElementMethodData,
                       CoupledNeumannEntity,
                       FiniteElementModule>
GalerkinStructMechHeat3DDiffusiveDispCoupledNeumannEntity("GalerkinStructMechHeat3DDiffusiveDisp");

MethodStrategyProvider<GalerkinCoupledNeumannEntity,
                       FiniteElementMethodData,
                       CoupledNeumannEntity,
                       FiniteElementModule>
GalerkinStructMech3DDiffusiveDispCoupledNeumannEntity("GalerkinStructMech3DDiffusiveDisp");

MethodStrategyProvider<GalerkinCoupledNeumannEntity,
                       FiniteElementMethodData,
                       CoupledNeumannEntity,
                       FiniteElementModule>
GalerkinHeat3DDiffusivePrimCoupledNeumannEntity("GalerkinHeat3DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

GalerkinCoupledNeumannEntity::GalerkinCoupledNeumannEntity(const std::string& name) :
  CoupledNeumannEntity(name)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinCoupledNeumannEntity::~GalerkinCoupledNeumannEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinCoupledNeumannEntity::setup()
{
  CoupledNeumannEntity::setup();

}

//////////////////////////////////////////////////////////////////////////////

RealVector& GalerkinCoupledNeumannEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
//   const std::vector<Framework::State*>& vars = *(_localElemData->solValues);
  const RealVector& shapeF = (*_localElemData->solShapeF)[iQuadPoint];
//   const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];
  Framework::GeometricEntity* const geo = _localElemData->cell;

  RealVector& state = *((*_localElemData->solValues)[iQuadPoint]);
  RealVector& coord = *((*_localElemData->coord)[iQuadPoint]);

  cf_assert(_isAcceptedIdx < _isAccepted.size());
  if((_isAccepted[_isAcceptedIdx]) >= 0.)
  {
    CFout << "Accepted: " << _isAcceptedIdx << "/" << _isAccepted.size()  << "\n";
    CFout << "Data: " << _index << "/" << _interfaceData.size() << "\n";
    CFout << "Coord: " << coord << "\n";
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
    else
    {
      _result = _interfaceData[_index];
      _index++;
    }
  }
  else // not accepted
  {
    /// @note why not derive from NeumannEntity and if not accepted, use NeumannEntity
    CFout << "Refused: " << _isAcceptedIdx << "\n";
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
  _result *= shapeF[iState];
  _isAcceptedIdx++;
  return _result;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

