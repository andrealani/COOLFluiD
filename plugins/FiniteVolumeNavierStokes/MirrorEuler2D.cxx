#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/MirrorEuler2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MirrorEuler2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
mirrorEuler2DFVMCCProvider("MirrorEuler2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

MirrorEuler2D::MirrorEuler2D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

MirrorEuler2D::~MirrorEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void MirrorEuler2D::setup()
{
  FVMCC_BC::setup();
	
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void MirrorEuler2D::setGhostState(GeometricEntity *const face)
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
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal vn = _dataInnerState[EulerTerm::VX]*nx +
    _dataInnerState[EulerTerm::VY]*ny;

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::VX] = _dataInnerState[EulerTerm::VX] - 2.0*vn*nx;
  _dataGhostState[EulerTerm::VY] = _dataInnerState[EulerTerm::VY] - 2.0*vn*ny;
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::P] = _dataInnerState[EulerTerm::P];
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
       + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*
       _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
  
  _dataGhostState[EulerTerm::T] = _dataInnerState[EulerTerm::T];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
