#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DProjectionPhotosphereGroth.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<MirrorMHD3DProjectionPhotosphereGroth, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DProjectionPhotosphereGrothFVMCCProvider("MirrorMHD3DProjectionPhotosphereGrothFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DProjectionPhotosphereGroth::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("rhoFixed", "Non-dimensional density value that is to be fixed in the ghost cells.");
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionPhotosphereGroth::MirrorMHD3DProjectionPhotosphereGroth(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  _cartesianSphericalTMInnerState(),
  _sphericalCartesianTMInnerState(),
  _cartesianSphericalTMGhostState(),
  _sphericalCartesianTMGhostState()
{
  addConfigOptionsTo(this);

  _rhoFixed = 1.0;
  setParameter("rhoFixed",&_rhoFixed);

}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionPhotosphereGroth::~MirrorMHD3DProjectionPhotosphereGroth()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPhotosphereGroth::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPhotosphereGroth::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _cartesianSphericalTMInnerState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
             Framework::PhysicalModelStack::getActive()->getDim());
  _sphericalCartesianTMInnerState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
                Framework::PhysicalModelStack::getActive()->getDim());
  _cartesianSphericalTMGhostState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
             Framework::PhysicalModelStack::getActive()->getDim());
  _sphericalCartesianTMGhostState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
                Framework::PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPhotosphereGroth::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

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

  RealVector innerStateCoordsSpherical(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector ghostStateCoordsSpherical(Framework::PhysicalModelStack::getActive()->getDim());

  const RealVector innerStateCoords = innerState->getCoordinates();
  const RealVector ghostStateCoords = ghostState->getCoordinates();

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  // set the transformation matrices between Cartesian and spherical coordinate systems
  _varSet->setTransformationMatrices(innerStateCoords,innerStateCoordsSpherical,_cartesianSphericalTMInnerState,_sphericalCartesianTMInnerState);
  _varSet->setTransformationMatrices(ghostStateCoords,ghostStateCoordsSpherical,_cartesianSphericalTMGhostState,_sphericalCartesianTMGhostState);

  RealVector VCartesianInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VSphericalInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VCartesianGhostState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VSphericalGhostState(Framework::PhysicalModelStack::getActive()->getDim());

  VCartesianInnerState[0] = _dataInnerState[MHDProjectionTerm::VX];
  VCartesianInnerState[1] = _dataInnerState[MHDProjectionTerm::VY];
  VCartesianInnerState[2] = _dataInnerState[MHDProjectionTerm::VZ];

  const CFreal vn = VCartesianInnerState[0]*nx + VCartesianInnerState[1]*ny + VCartesianInnerState[2]*nz;

  VSphericalInnerState = _cartesianSphericalTMInnerState*VCartesianInnerState;

  if (vn < 0.0) {
	// inflow through the boundary
	VSphericalGhostState[0] = innerStateCoordsSpherical[0]*innerStateCoordsSpherical[0]*_dataInnerState[MHDProjectionTerm::RHO]*VSphericalInnerState[0]/
			(ghostStateCoordsSpherical[0]*ghostStateCoordsSpherical[0]*_dataGhostState[MHDProjectionTerm::RHO]);
	VSphericalGhostState[1] = 0.0;
	VSphericalGhostState[2] = 0.0;
	VCartesianGhostState = _sphericalCartesianTMGhostState*VSphericalGhostState;	

        _dataGhostState[MHDProjectionTerm::RHO] = 2.*_rhoFixed - _dataInnerState[MHDProjectionTerm::RHO];
	_dataGhostState[MHDProjectionTerm::VX] = VCartesianGhostState[0];
        _dataGhostState[MHDProjectionTerm::VY] = VCartesianGhostState[1];
        _dataGhostState[MHDProjectionTerm::VZ] = VCartesianGhostState[2];
	_dataGhostState[MHDProjectionTerm::BX] = -_dataInnerState[MHDProjectionTerm::BX];
        _dataGhostState[MHDProjectionTerm::BY] = -_dataInnerState[MHDProjectionTerm::BY];
        _dataGhostState[MHDProjectionTerm::BZ] = -_dataInnerState[MHDProjectionTerm::BZ];
        _dataGhostState[MHDProjectionTerm::V] = fabs(VSphericalGhostState[0]);
	_dataGhostState[MHDProjectionTerm::P] = 1.0/_varSet->getModel()->getGamma();
        _dataGhostState[MHDProjectionTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDProjectionTerm::P]/_dataGhostState[MHDProjectionTerm::RHO]);
        _dataGhostState[MHDProjectionTerm::B] = _dataInnerState[MHDProjectionTerm::B];
        _dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];
  }
  if (vn >= 0.0) {
	// outflow through the boundary
	_dataGhostState[MHDProjectionTerm::RHO] = _dataInnerState[MHDProjectionTerm::RHO];
	_dataGhostState[MHDProjectionTerm::VX] = -_dataInnerState[MHDProjectionTerm::VX];
        _dataGhostState[MHDProjectionTerm::VY] = -_dataInnerState[MHDProjectionTerm::VY];
        _dataGhostState[MHDProjectionTerm::VZ] = -_dataInnerState[MHDProjectionTerm::VZ];
        _dataGhostState[MHDProjectionTerm::BX] = -_dataInnerState[MHDProjectionTerm::BX];
        _dataGhostState[MHDProjectionTerm::BY] = -_dataInnerState[MHDProjectionTerm::BY];
        _dataGhostState[MHDProjectionTerm::BZ] = -_dataInnerState[MHDProjectionTerm::BZ];
        _dataGhostState[MHDProjectionTerm::V] = _dataInnerState[MHDProjectionTerm::V];
	_dataGhostState[MHDProjectionTerm::P] = _dataInnerState[MHDProjectionTerm::P];
	_dataGhostState[MHDProjectionTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDProjectionTerm::P]/_dataGhostState[MHDProjectionTerm::RHO]);
	_dataGhostState[MHDProjectionTerm::B] = _dataInnerState[MHDProjectionTerm::B];
	_dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];
  }	
	
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
