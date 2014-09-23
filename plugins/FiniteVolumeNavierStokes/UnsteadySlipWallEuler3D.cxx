#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/UnsteadySlipWallEuler3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UnsteadySlipWallEuler3D, CellCenterFVMData, FiniteVolumeNavierStokesModule> UnsteadySlipWallEuler3DFVMCCProvider("UnsteadySlipWallEuler3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

UnsteadySlipWallEuler3D::UnsteadySlipWallEuler3D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes"),
  _bCoord(),
  _speed()
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySlipWallEuler3D::~UnsteadySlipWallEuler3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySlipWallEuler3D::setup()
{

  FVMCC_BC::setup();

  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  _speed.resize(PhysicalModelStack::getActive()->getDim());

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySlipWallEuler3D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = normals[startID + 2];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  computeGhostStateSpeed(face);

  const CFreal vn = _dataInnerState[EulerTerm::VX]*nx +
                    _dataInnerState[EulerTerm::VY]*ny+
                    _dataInnerState[EulerTerm::VZ]*nz;

  const CFreal normalSpeed = _speed[XX]*nx + _speed[YY]*ny + _speed[ZZ]*nz;

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.);

  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::VX] = _dataInnerState[EulerTerm::VX] - 2.0*vn*nx + 2.0*normalSpeed*nx;
  _dataGhostState[EulerTerm::VY] = _dataInnerState[EulerTerm::VY] - 2.0*vn*ny + 2.0*normalSpeed*ny;
  _dataGhostState[EulerTerm::VZ] = _dataInnerState[EulerTerm::VZ] - 2.0*vn*nz + 2.0*normalSpeed*nz;
  
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				       _dataGhostState[EulerTerm::VX] +
				       _dataGhostState[EulerTerm::VY]*
				       _dataGhostState[EulerTerm::VY] +
				       _dataGhostState[EulerTerm::VZ]*
				       _dataGhostState[EulerTerm::VZ]);
  
  ///@todo modify the pressure to account for acceleration of the flow!!!!!!!!!
  // dp/dn = -rho*n*a_w

  const CFreal pressureCorrection = 0.;
  _dataGhostState[EulerTerm::P] = _dataInnerState[EulerTerm::P] + pressureCorrection;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
       + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*
       _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = _dataInnerState[EulerTerm::T];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadySlipWallEuler3D::computeGhostStateSpeed(GeometricEntity *const face)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFuint nbNodes = face->nbNodes();

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();

  ///Compute the speed of the projected point from the boundary nodes
  // The point N now lies on the face.
  // We have to interpolate the speed at N from the speed of the nodes.
  // Do an averaging weighted by the distance...
  RealVector sumSpeedOverDistance(dim);
  RealVector vector(dim);
  RealVector speedNode(dim);

  sumSpeedOverDistance = 0.;
  CFreal sumOverDistance = 0.;

  for(CFuint iNode=0;iNode < nbNodes;iNode++)
  {
    const CFuint nodeID = face->getNode(iNode)->getLocalID();
    vector = ((*(face->getNode(iNode))) - _bCoord);
    CFreal norm = vector.norm2();
    speedNode = (*(futureNodes[nodeID]) - *(pastNodes[nodeID]))/dt;

    for(CFuint iDim=0;iDim < dim;iDim++)
    {
      sumSpeedOverDistance[iDim] += speedNode[iDim]/(norm);
    }
    sumOverDistance += 1./norm;
  }

  for(CFuint iDim=0;iDim < dim;iDim++)
  {
    _speed[iDim] = sumSpeedOverDistance[iDim] / sumOverDistance;
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
UnsteadySlipWallEuler3D::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_BC::needsSockets();

  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
