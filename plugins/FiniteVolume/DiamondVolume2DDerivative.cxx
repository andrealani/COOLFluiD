#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "DiamondVolume2DDerivative.hh"
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

MethodStrategyProvider<DiamondVolume2DDerivative,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
diamondVolume2DDerivativeProvider("DiamondVolume2D");

//////////////////////////////////////////////////////////////////////////////

DiamondVolume2DDerivative::DiamondVolume2DDerivative(const std::string& name) :
  DerivativeComputer(name),
  _n01(),
  _n12(),
  _n23(),
  _n30(),
  _nodes(4),
  _midNode(),
  _weights(4)
{
}

//////////////////////////////////////////////////////////////////////////////

DiamondVolume2DDerivative::~DiamondVolume2DDerivative()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume2DDerivative::setup()
{
  DerivativeComputer::setup();

  _n01.resize(PhysicalModelStack::getActive()->getDim());
  _n12.resize(PhysicalModelStack::getActive()->getDim());
  _n23.resize(PhysicalModelStack::getActive()->getDim());
  _n30.resize(PhysicalModelStack::getActive()->getDim());
  _midNode.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume2DDerivative::computeGradients(const RealMatrix& values,
						 vector<RealVector*>& gradients)
{
  // Green Gauss is applied in the diamond volume to compute the
  // gradients
  const CFreal invVolume = 0.5/_volume;
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);

  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];

    const CFreal v01 = values(i,0) + values(i,1);
    const CFreal v12 = values(i,1) + values(i,2);
    const CFreal v23 = values(i,2) + values(i,3);
    const CFreal v30 = values(i,3) + values(i,0);

    for (CFuint ix = 0; ix < DIM_2D; ++ix) {
      grad[ix] = invVolume*(_n01[ix]*v01 + _n12[ix]*v12 +
			    _n23[ix]*v23 + _n30[ix]*v30);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume2DDerivative::computeControlVolume
(vector<RealVector*>& states, GeometricEntity *const geo)
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  states[0] = &nodalStates[geo->getNode(0)->getLocalID()];
  states[1] = geo->getState(0);
  states[2] = &nodalStates[geo->getNode(1)->getLocalID()];
  states[3] = geo->getState(1);
  
  // construct the quadrilateral volume and the normals around the face
  Node& node0 = *geo->getNode(0);
  Node& node1 = geo->getState(0)->getCoordinates();
  Node& node2 = *geo->getNode(1);
  Node& node3 = geo->getState(1)->getCoordinates();

  _n01[XX] = -node0[YY] + node1[YY];
  _n01[YY] = -node1[XX] + node0[XX];

  _n12[XX] = -node1[YY] + node2[YY];
  _n12[YY] = -node2[XX] + node1[XX];

  _n23[XX] = -node2[YY] + node3[YY];
  _n23[YY] = -node3[XX] + node2[XX];

  _n30[XX] = -node3[YY] + node0[YY];
  _n30[YY] = -node0[XX] + node3[XX];

  const CFreal det = (node2[XX] - node0[XX])*(node3[YY] - node1[YY]) -
    (node3[XX] - node1[XX])*(node2[YY] - node0[YY]);

  if (det < 0.0) {
    _n01 *= -1.;
    _n12 *= -1.;
    _n23 *= -1.;
    _n30 *= -1.;
  }

  // store the nodes of the control volume
  _nodes[0] = &node0;
  _nodes[1] = &node1;
  _nodes[2] = &node2;
  _nodes[3] = &node3;

  cf_assert(std::abs(det) > 0);
  _volume = 0.5*std::abs(det);
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume2DDerivative::computeAverageValues
(GeometricEntity *const geo,
 const vector<RealVector*>& values,
 RealVector& avValues)
{
  const CFuint nbValues = avValues.size();
  for (CFuint i = 0; i < nbValues; ++i) {
    avValues[i] = 0.5*((*values[0])[i] + (*values[2])[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<vector<RealVector> > DiamondVolume2DDerivative::getGradientsJacob()
{
  _gradientsJacob[0] = (0.5/_volume)*(_n01 + _n12);
  _gradientsJacob[1] = (0.5/_volume)*(_n30 + _n23);
  
  return &_gradientsJacob;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
