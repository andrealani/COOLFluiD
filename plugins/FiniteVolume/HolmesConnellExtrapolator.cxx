#include "MathTools/MatrixInverterT.hh"
#include "MathTools/MathConsts.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"

#include "FiniteVolume/HolmesConnellExtrapolator.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HolmesConnellExtrapolator,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeModule>
holmesConnellExtrapolatorProvider("HolmesConnell");

//////////////////////////////////////////////////////////////////////////////

HolmesConnellExtrapolator::HolmesConnellExtrapolator(const std::string& name) :
  NodalStatesExtrapolator<CellCenterFVMData>(name),
  _coeffs()
{
}

//////////////////////////////////////////////////////////////////////////////

HolmesConnellExtrapolator::~HolmesConnellExtrapolator()
{
}

//////////////////////////////////////////////////////////////////////////////

void HolmesConnellExtrapolator::setup()
{
  NodalStatesExtrapolator<CellCenterFVMData>::setup();
  
  // compute the reconstruction coefficients
  // (they only depend on the given mesh)
  computeCoeffs();
}

//////////////////////////////////////////////////////////////////////////////

void HolmesConnellExtrapolator::computeCoeffs()
{
  MathTools::MatrixInverterT<2> inverter2;
  MathTools::MatrixInverterT<3> inverter3;

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(dim == 2 || dim == 3);

  _coeffs.resize(dim*nbNodes);

  RealMatrix matrix(dim,dim);
  RealMatrix inverse(dim,dim);
  RealVector coeff(dim);
  RealVector rhs(dim);

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    // reset to 0 the matrix
    matrix = 0.0;
    coeff = 0.0;
    rhs = 0.0;

    const CFuint nbNeighborStates = _neighborStates[iNode].size();
    const Node& currNode = *nodes[iNode];

    if (dim == 3) {
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
	const State *const neighState = _neighborStates[iNode][iState];
	const Node& neighNode = neighState->getCoordinates();
	
	const CFreal dx = neighNode[0] - currNode[0];
	const CFreal dy = neighNode[1] - currNode[1];
	const CFreal dz = neighNode[2] - currNode[2];
	
	// update the matrix
	matrix(0,0) += dx*dx;
	matrix(0,1) += dx*dy;
	matrix(0,2) += dx*dz;
	
	matrix(1,0) += dy*dx;
	matrix(1,1) += dy*dy;
	matrix(1,2) += dy*dz;
	
	matrix(2,0) += dz*dx;
	matrix(2,1) += dz*dy;
	matrix(2,2) += dz*dz;
	
	// update the rhs
	rhs[0] -= dx;
	rhs[1] -= dy;
	rhs[2] -= dz;
      }
      
      inverter3.invert(matrix,inverse);
    }
    else {
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
	const State *const neighState = _neighborStates[iNode][iState];
	const Node& neighNode = neighState->getCoordinates();
	
	const CFreal dx = neighNode[0] - currNode[0];
	const CFreal dy = neighNode[1] - currNode[1];
	
	// update the matrix
	matrix(0,0) += dx*dx;
	matrix(0,1) += dx*dy;
	
	matrix(1,0) += dy*dx;
	matrix(1,1) += dy*dy;
	
	// update the rhs
	rhs[0] -= dx;
	rhs[1] -= dy;
      }
      
      inverter2.invert(matrix,inverse);
    }
    
    // compute the coeffs
    coeff = inverse*rhs;
    
    // set the coeffs
    const CFuint startID = iNode*dim;
    for (CFuint i = 0; i < dim; ++i) {
      _coeffs[startID + i] = coeff[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void HolmesConnellExtrapolator::extrapolateInAllNodes()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(dim == 2 || dim == 3);

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const CFuint nbNeighborStates = _neighborStates[iNode].size();
    const Node& currNode = *nodes[iNode];
    const CFuint startID = iNode*dim;

    // reset to 0 the nodal states and weights
    nodalStates[iNode] = 0.;
    CFreal sumWeights = 0.0;

    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const State *const neighState = _neighborStates[iNode][iState];
      const Node& neighNode = neighState->getCoordinates();

      CFreal weight = 1.0;
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
  const CFreal dx = neighNode[iDim] - currNode[iDim];
  weight += _coeffs[startID + iDim]*dx;
      }

      nodalStates[iNode] += (*neighState)*weight;
      sumWeights += weight;
    }

    const bool isValid = (std::abs(sumWeights) > MathTools::MathConsts::CFrealEps());
    if (isValid) {
      nodalStates[iNode] *= 1./sumWeights;
    }
    else {
      // if the sum of weights is too small (< eps) use ditsnace based
      // reconstruction which is never singular
      nodalStates[iNode] = 0.0;

      CFreal sumDr = 0.0;
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
  const State *const neighState = _neighborStates[iNode][iState];
  const Node& neighNode = neighState->getCoordinates();
  const CFreal dr =  1./MathTools::MathFunctions::getDistance(neighNode, currNode);
  nodalStates[iNode] += (*neighState)*dr;
  sumDr += dr;
      }

      nodalStates[iNode] *= 1./sumDr;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void HolmesConnellExtrapolator::extrapolateInNodes
(const vector<Node*>& nodes)
{

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

  const CFuint nbNodes = nodes.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const Node& currNode = *nodes[iNode];
    const CFuint nodeID = currNode.getLocalID();
    const CFuint nbNeighborStates = _neighborStates[nodeID].size();
    const CFuint startID = nodeID*dim;

    // reset to 0 the nodal states
    nodalStates[nodeID] = 0.;
    CFreal sumWeights = 0.0;

    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const State *const neighState = _neighborStates[nodeID][iState];
      const Node& neighNode = neighState->getCoordinates();

      CFreal weight = 1.0;
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
  const CFreal dx = neighNode[iDim] - currNode[iDim];
  weight += _coeffs[startID + iDim]*dx;
      }

      nodalStates[nodeID] += (*neighState)*weight;
      sumWeights += weight;
    }

    const bool isValid = (std::abs(sumWeights) > MathTools::MathConsts::CFrealEps());
    if (isValid) {
      nodalStates[nodeID] *= 1./sumWeights;
    }
    else {
      // if the sum of weights is too small (< eps) use ditsnace based
      // reconstruction which is never singular
      nodalStates[nodeID] = 0.0;

      CFreal sumDr = 0.0;
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
  const State *const neighState = _neighborStates[nodeID][iState];
  const Node& neighNode = neighState->getCoordinates();
  const CFreal dr =  1./MathTools::MathFunctions::getDistance(neighNode, currNode);
  nodalStates[nodeID] += (*neighState)*dr;
  sumDr += dr;
      }

      nodalStates[nodeID] *= 1./sumDr;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
