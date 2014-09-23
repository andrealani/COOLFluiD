#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNavierStokes/MirrorVelocity3DTurb.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

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

MethodCommandProvider<MirrorVelocity3DTurb, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
mirrorVelocity3DTurbFVMCCProvider("MirrorVelocity3DTurbFVMCC");

//////////////////////////////////////////////////////////////////////////////

MirrorVelocity3DTurb::MirrorVelocity3DTurb(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

MirrorVelocity3DTurb::~MirrorVelocity3DTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void MirrorVelocity3DTurb::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = normals[startID + 2];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  const CFreal vn = _dataInnerState[EulerTerm::VX]*nx +
                    _dataInnerState[EulerTerm::VY]*ny +
                    _dataInnerState[EulerTerm::VZ]*nz;

  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;

  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::VX]  = _dataInnerState[EulerTerm::VX] - 2.0*vn*nx;
  _dataGhostState[EulerTerm::VY]  = _dataInnerState[EulerTerm::VY] - 2.0*vn*ny;
  _dataGhostState[EulerTerm::VZ]  = _dataInnerState[EulerTerm::VZ] - 2.0*vn*nz;
  _dataGhostState[EulerTerm::V]   = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::P]   = _dataInnerState[EulerTerm::P];
  _dataGhostState[EulerTerm::H]   = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
                                    + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*
                                          _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
                                             _dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::T] = _dataInnerState[EulerTerm::T];
  
  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);

  for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
  {
    _dataGhostState[iK + iTurb] = _dataInnerState[iK + iTurb];
  }

  // set the ghost state starting from the physical data
  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void MirrorVelocity3DTurb::setup()
{
  FVMCC_BC::setup();
 
  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();
  cf_assert(_varSetTurb.isNotNull());

  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  if(_turbVars.size() == 0){
    _turbVars.resize(_varSetTurb->getModel()->getNbScalarVars(0));
  }
  
  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb3DVarSet>();
  cf_assert(_diffVarTurb.isNotNull());
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
