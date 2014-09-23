#include "StateDiffDerivative.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StateDiffDerivative,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
stateDiffDerivativeProvider("StateDiff");

//////////////////////////////////////////////////////////////////////////////

StateDiffDerivative::StateDiffDerivative(const std::string& name) :
  DerivativeComputer(name),
  socket_volumes("volumes"),
  _dr(0.0),
  _vectorLR()
{
}

//////////////////////////////////////////////////////////////////////////////

StateDiffDerivative::~StateDiffDerivative()
{
}

//////////////////////////////////////////////////////////////////////////////

void StateDiffDerivative::setup()
{
  DerivativeComputer::setup();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _vectorLR.resize(dim);
}

//////////////////////////////////////////////////////////////////////////////

void StateDiffDerivative::computeGradients(const RealMatrix& values,
					   vector<RealVector*>& gradients)
{
  // Green Gauss is applied in the diamond volume to compute the
  // gradients
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    const CFreal dv = (values(i,3) - values(i,1))/_dr;
    for (CFuint ix = 0; ix < dim; ++ix) {
      grad[ix] = dv*_vectorLR[ix]; // + o - according to the outward
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StateDiffDerivative::computeControlVolume
(std::vector<RealVector*>& states, GeometricEntity *const geo)
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  states[0] = &nodalStates[geo->getNode(0)->getLocalID()];
  states[1] = geo->getState(0);
  states[2] = &nodalStates[geo->getNode(1)->getLocalID()];
  states[3] = geo->getState(1);
  
  //   Node& node0 = *geo->getNode(0);
  Node& node1 = geo->getState(0)->getCoordinates();
  //   Node& node2 = *geo->getNode(1);
  Node& node3 = geo->getState(1)->getCoordinates();
  _dr = MathTools::MathFunctions::getDistance(node1, node3);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint ix = 0; ix < dim; ++ix) {
    _vectorLR[ix] =  node3[ix] - node1[ix];
  }
  _vectorLR /= _dr;
  
  // handle to ID of the cell for which the normal is outward
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  const CFuint faceID = geo->getID();
  if (static_cast<CFuint>(isOutward[faceID]) != geo->getState(0)->getLocalID()) {
    _vectorLR *= -1.;
  }
  
  //  _volume = _dr*MathTools::MathFunctions::getDistance(node0, node2);
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  if (!geo->getState(1)->isGhost()) {
    _volume = 0.5*(volumes[geo->getState(0)->getLocalID()] + volumes[geo->getState(1)->getLocalID()]);
  }
  else {
    _volume = volumes[geo->getState(0)->getLocalID()];
  }
  
  cf_assert(_volume > 0.0);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void StateDiffDerivative::computeAverageValues
(GeometricEntity *const geo,const vector<RealVector*>& values,RealVector& avValues)
{
  const CFuint nbValues = avValues.size();
  for (CFuint i = 0; i < nbValues; ++i) {
    avValues[i] = 0.5*((*values[0])[i] + (*values[2])[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<vector<RealVector> > StateDiffDerivative::getGradientsJacob()
{
  /// @TODO to be implemented properly
  return &_gradientsJacob;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StateDiffDerivative::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    DerivativeComputer::needsSockets();

  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
