#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/UnsteadyNoSlipWallAdiabaticNSTurb2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

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

MethodCommandProvider<UnsteadyNoSlipWallAdiabaticNSTurb2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> UnsteadyNoSlipWallAdiabaticNSTurb2DFVMCCProvider("UnsteadyNoSlipWallAdiabaticNSTurb2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

UnsteadyNoSlipWallAdiabaticNSTurb2D::UnsteadyNoSlipWallAdiabaticNSTurb2D(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes"),
  _coord(),
  _speed(),
  socket_wallDistance("wallDistance")
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyNoSlipWallAdiabaticNSTurb2D::~UnsteadyNoSlipWallAdiabaticNSTurb2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNoSlipWallAdiabaticNSTurb2D::setup()
{

  FVMCC_BC::setup();

  _coord.resize(PhysicalModelStack::getActive()->getDim());
  _speed.resize(PhysicalModelStack::getActive()->getDim());
  
  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb2DVarSet>();
  
  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNoSlipWallAdiabaticNSTurb2D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  computeGhostStateSpeed(face);

//unused//  const CFreal vn = _dataInnerState[EulerTerm::VX]*nx + _dataInnerState[EulerTerm::VY]*ny;
  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal ghostT = _dataInnerState[EulerTerm::P]/
    (R*_dataInnerState[EulerTerm::RHO]);

  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.);

//   if(((*_futureNodes[nodeID0])[YY] < 0.06)&&((*_futureNodes[nodeID0])[XX] < 0.05)&&((*_futureNodes[nodeID0])[XX] > 0.)){
// std::cout << "------------------------------------" << std::endl;
// std::cout << " " << std::endl;
//   DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();
//   const CFuint nodeID0 = face->getNode(0)->getLocalID();
//   const CFuint nodeID1 = face->getNode(1)->getLocalID();

//std::cout << "Face Nodes: " << (*futureNodes[nodeID0]) << std::endl;
//std::cout << "            " << (*futureNodes[nodeID1]) << std::endl;
//   //std::cout << "GhostNode: " << ghostState->getCoordinates() << std::endl;
//   //std::cout << "InnerNode: " << innerState->getCoordinates() << std::endl;
//  std::cout << "Speed: " << speed << std::endl;
//   std::cout << "Normal: " << nx << " " << ny << std::endl;
//std::cout << "normalSpeed: " << normalSpeed << std::endl;
//   }

  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::VX]  = -_dataInnerState[EulerTerm::VX] + 2.0*_speed[XX];
  _dataGhostState[EulerTerm::VY]  = -_dataInnerState[EulerTerm::VY] + 2.0*_speed[YY];

  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				       _dataGhostState[EulerTerm::VX] +
				       _dataGhostState[EulerTerm::VY]*
				       _dataGhostState[EulerTerm::VY]);

  ///@todo modify the pressure to account for acceleration of the flow!!!!!!!!!
  // dp/dn = -rho*n*a_w

  const CFreal pressureCorrection = 0.;
  _dataGhostState[EulerTerm::P] = _dataInnerState[EulerTerm::P] + pressureCorrection;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
       + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*
       _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = ghostT;
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] -
     (_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);



  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);

  //k=0 at the wall
  //(in reality, put a very small value to avoid negative values close to the wall)
  const CFreal kWall = 1.e-20;
  _dataGhostState[iK] = kWall - _dataInnerState[iK];

  //At the wall, according to Menter's definition, Omega_w = 10*((6*nu)/(beta1 * y0 *y0))
  if(nbTurbVars == 2){

      ///Compute y0
      DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
      CFreal y0 = wallDistance[innerState->getLocalID()];

      ///This is a temporary fix for the flat plate because, at initialization, the wall Distance is not known!!!!
      ///@todo remove this line...
      if(SubSystemStatusStack::getActive()->getNbIter() == 0) y0 = innerState->getCoordinates()[YY];

      //avoid too small distances
      y0 = max(y0, 10.e-10);

      CFreal pdim = _dataInnerState[EulerTerm::P] * _varSetTurb->getModel()->getPressRef();
      CFreal Tdim = _dataInnerState[EulerTerm::T] * _varSetTurb->getModel()->getTempRef();
      
      CFreal mu = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
	_diffVarTurb->getModel().getReferencePhysicalData()[NSTurbTerm::MU];
      
      CFreal nu = mu / _dataInnerState[EulerTerm::RHO];

      //this is not the best, but it avoids having to code another BC! because I
      //would have to dynamic cast to the KOmega varset to get the beta1
      CFreal beta1 = 0.075;

      ///@todo here should this be adimensionalized (by the distance)???
      //Menter's definition
      CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);

      //Wilcox's definition
      //CFreal omegaWall = (1. * 6. * nu) / (beta1 * y0 * y0);

      _dataGhostState[iK + 1] = 2.0*omegaWall - _dataInnerState[iK + 1];
  }


  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNoSlipWallAdiabaticNSTurb2D::computeGhostStateSpeed(GeometricEntity *const face)
{

  ///Compute the projection of the cell center on the boundary face
  // The equation of the plane containing the boundary face
  // and the given node (xp, yp, zp) (first node of the face),
  // with normal (a,b,c) is
  // a*x + b*y + c*z + k = 0
  // with k = -a*xp - b*yp - c*zp
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
//unused//  const CFuint nbNodes = face->nbNodes();
  const Node* const node = face->getNode(0);
  const State* const state = face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*dim;
  RealVector faceNormal(dim);

  DataHandle< CFreal> normals = socket_normals.getDataHandle();

  // set the current normal
  for (CFuint iDim = 0; iDim < dim; ++iDim) {
    faceNormal[iDim] = normals[startID + iDim];
  }

  CFLogDebugMax( "faceNormal = " << faceNormal << "\n");
  CFLogDebugMax( "state coord = " << state->getCoordinates() << "\n");

  const CFreal k = - MathFunctions::innerProd(faceNormal, *node);

  // t is parameter for vectorial representation of a straight line
  // in space
  // t = (a*xM + b*yM + c*zM + k)/(a*a + b*b + c*c)

  const CFreal n2 = MathFunctions::innerProd(faceNormal, faceNormal);

  cf_assert(std::abs(n2) > 0.0);

  const RealVector& stateCoord = state->getCoordinates();
  const CFreal t = (MathFunctions::innerProd(faceNormal,stateCoord) + k)/n2;

  // The point N, projection from M (position of the other state,
  // internal neighbor of the face) with respect to the
  // given plane is given by
  // (xN, yN, zN) = (xM, yM, zM) - t*(a,b,c)

  for (CFuint iDim = 0; iDim < dim; ++iDim) {
    _coord[iDim] =  stateCoord[iDim] - t*faceNormal[iDim];
  }

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();

  ///Compute the speed of the projected point from the boundary nodes
  // The point N now lies on the face.
  // We have to interpolate the speed at N from the speed of the nodes.
  // Do an averaging weighted by the distance...

  // for 2D just compute more easily
  const CFuint nodeID0 = face->getNode(0)->getLocalID();
  const CFuint nodeID1 = face->getNode(1)->getLocalID();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  const RealVector speedNode0 = (*(futureNodes[nodeID0]) - *(pastNodes[nodeID0]))/dt;
  const RealVector speedNode1 = (*(futureNodes[nodeID1]) - *(pastNodes[nodeID1]))/dt;

  const RealVector projectionToNode0 = ((*(face->getNode(0))) - _coord);
  const RealVector node1ToNode0 = ((*(face->getNode(1))) - (*(face->getNode(0))));
  const CFreal ratio = projectionToNode0.norm2()/node1ToNode0.norm2();
// std::cout << "Coord: " << _coord << std::endl;
// std::cout << "Node0: " << (*(face->getNode(0))) << std::endl;
// std::cout << "Node1: " << (*(face->getNode(1))) << std::endl;

  _speed = (ratio * (speedNode1 - speedNode0));
  _speed += speedNode0;

/*  RealVector sumCoordOverDistance(dim);
  RealVector sumOverDistance(dim);

  sumCoordOverDistance = 0.;
  sumOverDistance = 0.;

  for(CFuint iNode=0;iNode < nbNodes;iNode++)
  {
    for(CFuint iDim=0;iDim < dim;iDim++)
    {
      sumCoordOverDistance[iDim] += (*(face->getNode(iNode)))[iDim] / ((*(face->getNode(iNode)))[iDim] - _coord[iDim]);
      sumOverDistance[iDim] += 1./((*(face->getNode(iNode)))[iDim] - _coord[iDim]);
    }
  }

  for(CFuint iDim=0;iDim < dim;iDim++)
  {
    _speed[iDim] = sumCoordOverDistance[iDim] / sumOverDistance[iDim];
  }*/
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
UnsteadyNoSlipWallAdiabaticNSTurb2D::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_BC::needsSockets();

  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
