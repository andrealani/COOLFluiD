#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DProjectionPolytropicPhotosphereGroth.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionPolytropicVarSet.hh"
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

MethodCommandProvider<MirrorMHD3DProjectionPolytropicPhotosphereGroth, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DProjectionPolytropicPhotosphereGrothFVMCCProvider("MirrorMHD3DProjectionPolytropicPhotosphereGrothFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DProjectionPolytropicPhotosphereGroth::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("rhoFixed", "Non-dimensional density value that is to be fixed in the ghost cells.");
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionPolytropicPhotosphereGroth::MirrorMHD3DProjectionPolytropicPhotosphereGroth(const std::string& name) :
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

MirrorMHD3DProjectionPolytropicPhotosphereGroth::~MirrorMHD3DProjectionPolytropicPhotosphereGroth()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPolytropicPhotosphereGroth::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPolytropicPhotosphereGroth::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionPolytropicVarSet>();
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

void MirrorMHD3DProjectionPolytropicPhotosphereGroth::setGhostState(GeometricEntity *const face)
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

  VCartesianInnerState[0] = _dataInnerState[MHDProjectionPolytropicTerm::VX];
  VCartesianInnerState[1] = _dataInnerState[MHDProjectionPolytropicTerm::VY];
  VCartesianInnerState[2] = _dataInnerState[MHDProjectionPolytropicTerm::VZ];

  const CFreal vn = VCartesianInnerState[0]*nx + VCartesianInnerState[1]*ny + VCartesianInnerState[2]*nz;

  VSphericalInnerState = _cartesianSphericalTMInnerState*VCartesianInnerState;

  if (vn < 0.0) {
	// inflow through the boundary
	VSphericalGhostState[0] = innerStateCoordsSpherical[0]*innerStateCoordsSpherical[0]*_dataInnerState[MHDProjectionPolytropicTerm::RHO]*VSphericalInnerState[0]/
			(ghostStateCoordsSpherical[0]*ghostStateCoordsSpherical[0]*_dataGhostState[MHDProjectionPolytropicTerm::RHO]);
	VSphericalGhostState[1] = 0.0;
	VSphericalGhostState[2] = 0.0;
	VCartesianGhostState = _sphericalCartesianTMGhostState*VSphericalGhostState;	

        _dataGhostState[MHDProjectionPolytropicTerm::RHO] = 2.*_rhoFixed - _dataInnerState[MHDProjectionPolytropicTerm::RHO];
	_dataGhostState[MHDProjectionPolytropicTerm::VX] = VCartesianGhostState[0];
        _dataGhostState[MHDProjectionPolytropicTerm::VY] = VCartesianGhostState[1];
        _dataGhostState[MHDProjectionPolytropicTerm::VZ] = VCartesianGhostState[2];
	_dataGhostState[MHDProjectionPolytropicTerm::BX] = -_dataInnerState[MHDProjectionPolytropicTerm::BX];
        _dataGhostState[MHDProjectionPolytropicTerm::BY] = -_dataInnerState[MHDProjectionPolytropicTerm::BY];
        _dataGhostState[MHDProjectionPolytropicTerm::BZ] = -_dataInnerState[MHDProjectionPolytropicTerm::BZ];
        _dataGhostState[MHDProjectionPolytropicTerm::V] = fabs(VSphericalGhostState[0]);
	_dataGhostState[MHDProjectionPolytropicTerm::P] = pow(_dataGhostState[MHDProjectionPolytropicTerm::RHO],_dataGhostState[MHDProjectionPolytropicTerm::N])/_dataGhostState[MHDProjectionPolytropicTerm::N]; 
        _dataGhostState[MHDProjectionPolytropicTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDProjectionPolytropicTerm::P]/_dataGhostState[MHDProjectionPolytropicTerm::RHO]);
        _dataGhostState[MHDProjectionPolytropicTerm::B] = _dataInnerState[MHDProjectionPolytropicTerm::B];
        _dataGhostState[MHDProjectionPolytropicTerm::PHI] = _dataInnerState[MHDProjectionPolytropicTerm::PHI];
  }
  if (vn >= 0.0) {
	// outflow through the boundary
	_dataGhostState[MHDProjectionPolytropicTerm::RHO] = _dataInnerState[MHDProjectionPolytropicTerm::RHO];
	_dataGhostState[MHDProjectionPolytropicTerm::VX] = -_dataInnerState[MHDProjectionPolytropicTerm::VX];
        _dataGhostState[MHDProjectionPolytropicTerm::VY] = -_dataInnerState[MHDProjectionPolytropicTerm::VY];
        _dataGhostState[MHDProjectionPolytropicTerm::VZ] = -_dataInnerState[MHDProjectionPolytropicTerm::VZ];
        _dataGhostState[MHDProjectionPolytropicTerm::BX] = -_dataInnerState[MHDProjectionPolytropicTerm::BX];
        _dataGhostState[MHDProjectionPolytropicTerm::BY] = -_dataInnerState[MHDProjectionPolytropicTerm::BY];
        _dataGhostState[MHDProjectionPolytropicTerm::BZ] = -_dataInnerState[MHDProjectionPolytropicTerm::BZ];
        _dataGhostState[MHDProjectionPolytropicTerm::V] = _dataInnerState[MHDProjectionPolytropicTerm::V];
	_dataGhostState[MHDProjectionPolytropicTerm::P] = _dataInnerState[MHDProjectionPolytropicTerm::P];
	_dataGhostState[MHDProjectionPolytropicTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDProjectionPolytropicTerm::P]/_dataGhostState[MHDProjectionPolytropicTerm::RHO]);
	_dataGhostState[MHDProjectionPolytropicTerm::B] = _dataInnerState[MHDProjectionPolytropicTerm::B];
	_dataGhostState[MHDProjectionPolytropicTerm::PHI] = _dataInnerState[MHDProjectionPolytropicTerm::PHI];
  }	
	
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
