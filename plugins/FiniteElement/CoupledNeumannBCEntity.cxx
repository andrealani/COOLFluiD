#include "CoupledNeumannBCEntity.hh"
#include "Framework/IntegrableEntity.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

CoupledNeumannBCEntity::CoupledNeumannBCEntity() :
IntegrableEntity(),
_result(),
_index(0),
_isAcceptedIdx(0),
_iState(0),
_interfaceData(CFNULL),
_isAccepted(CFNULL),
_dim(),
_nbEq(),
_vars(),
_normal(),
_isRobin(false)
{
}

//////////////////////////////////////////////////////////////////////////////

CoupledNeumannBCEntity::~CoupledNeumannBCEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannBCEntity::setup()
{
  _result.resize(PhysicalModelStack::getActive()->getNbEq());
  _dim = PhysicalModelStack::getActive()->getDim();
  _nbEq = PhysicalModelStack::getActive()->getNbEq();
  _vars.resize(PhysicalModelStack::getActive()->getDim() + 1 + PhysicalModelStack::getActive()->getNbEq() + PhysicalModelStack::getActive()->getDim());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  _varsRobin.resize(2*PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

RealVector& CoupledNeumannBCEntity::operator() (const State& state, const RealVector& shapeF, GeometricEntity* geo)
{
  std::cout << FromHere().str() << std::endl;

  cf_assert(_isAcceptedIdx < _isAccepted.size());
  if((_isAccepted[_isAcceptedIdx]) >= 0.)
  {
/* CFout << "Data: " << _index << "/" << _interfaceData.size() << "\n";
 CFout << "Accepted: " << _isAcceptedIdx << "/" << _isAccepted.size()  << "\n";*/
 //CFout << "Coord: " << state.getCoordinates() << "\n";
    if(_isRobin)
    {

// Average value to avoid oscillations...
/*      std::vector<State*>* states = geo->getStates();
      const CFuint nbFaceStates = states->size();
      RealVector averageState = *((*states)[0]);
      for(CFuint i = 1; i < nbFaceStates; ++i) {
        averageState += *((*states)[i]);
      }
      averageState /= nbFaceStates;*/
      for(CFuint iEq = 0; iEq < _nbEq; ++iEq) {
        _varsRobin[iEq] = state[iEq];
      }
// CFout << "Data size: " << _interfaceData[_index].size() << "\n";
      for(CFuint iData = 0; iData < _interfaceData[_index].size(); ++iData) {
        _varsRobin[_nbEq + iData ] = (_interfaceData[_index])[iData];
      }
      std::cout << "Evaluating RobinBC \n";
      _vfRobin->evaluate(_varsRobin,_result);
      std::cout << "Coord: " << state.getCoordinates() << "\n";
      std::cout << "Temperature: " << state << "\n";
      std::cout << "Data: " << _interfaceData[_index] << "\n";
      std::cout << "Flux: " << _result << "\n";
      std::cout << std::endl;
      _index++;
    }
    else{
// CFout << "Coord: " << state.getCoordinates() << " - InterfaceData: " << _interfaceData[_index]  << "\n";
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
      throw Common::NotImplementedException("Robin BC is not implemented for not accepted points");
    }
    else
    {
//       _result = 0.;
//       CFout << "Not accepted : " << state.getCoordinates() << "\n";
      // Variables are Coord, Time, States, Normals
      // Set coordinates
      Node& coord = state.getCoordinates();
      for(CFuint i = 0; i < _dim; ++i) {
        _vars[i] = coord[i];
      }

      // Set time
      _vars[_dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

      // Set state variables
      for(CFuint i = 0; i < _nbEq; ++i) {
        _vars[i+_dim+1] = state[i];
      }

      // Set normals
      _normal = geo->computeAvgCellNormal();
      _normal.normalize();
      for(CFuint i = 0; i < _nbEq; ++i) {
        _vars[i + _nbEq + _dim + 1] = _normal[i];
      }

      _vf->evaluate(_vars,_result);
    }
  }
// CFout << "result * ShapeF[" << _iState << "]\n";
  _result *= shapeF[_iState];
  _isAcceptedIdx++;
  return _result;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace COOLFluiD

  } // namespace Numerics

} // namespace FiniteElement

//////////////////////////////////////////////////////////////////////////////

