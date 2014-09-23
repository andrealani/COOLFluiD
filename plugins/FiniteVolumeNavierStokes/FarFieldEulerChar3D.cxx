#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/FarFieldEulerChar3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerVarSet.hh"
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

MethodCommandProvider<FarFieldEulerChar3D, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
FarFieldEulerChar3DFVMCCProvider("FarFieldEulerChar3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

FarFieldEulerChar3D::FarFieldEulerChar3D(const std::string& name) :
  FarFieldEuler3D(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FarFieldEulerChar3D::~FarFieldEulerChar3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void FarFieldEulerChar3D::setGhostState(GeometricEntity *const face)
{
   State *const innerState = face->getState(0);
   State *const ghostState = face->getState(1);

   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.0;
   const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
   const CFreal R = _varSet->getModel()->getR();

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

   const CFreal u = _dataInnerState[EulerTerm::VX];
   const CFreal v = _dataInnerState[EulerTerm::VY];
   const CFreal w = _dataInnerState[EulerTerm::VZ];
   const CFreal rho = _dataInnerState[EulerTerm::RHO];
   const CFreal vn = u*nx + v*ny + w*nz;
   const CFreal pInnerState = _dataInnerState[EulerTerm::P];
   const CFreal aInnerState = _dataInnerState[EulerTerm::A];
   const CFreal machInner = vn / aInnerState;

   const CFreal vnInf = _uInf*nx + _vInf*ny + _wInf*nz;
   const CFreal rhoInf = _pressure/(R*_temperature);
   const CFreal aInf = sqrt(gamma*_pressure/rhoInf);
   //Compute the Riemann invariants
   const CFreal RPlus_n = vn + (2.*aInnerState)/(gammaMinus1);
   const CFreal RMin_n = vnInf  - (2.*aInf)/(gammaMinus1);

   CFreal RMinRmax = 0.5 * (RPlus_n+ RMin_n);

   // depending on the sign and magnitude of the local Mach number,
   // number of variables to be specified are determined

  // supersonic outlet case
  if (machInner >= 1.0)
  {
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = (*innerState)[1];
    (*ghostState)[2] = (*innerState)[2];
    (*ghostState)[3] = (*innerState)[3];
    (*ghostState)[4] = (*innerState)[4];
  }

  // supersonic inlet case
  if (machInner <= -1.0)
  {
     const CFreal rhoInf = _pressure/(R*_temperature);

     // set all the physical data corresponding to the ghost state
     _dataGhostState[EulerTerm::RHO] = 2.0*rhoInf - rho;
     _dataGhostState[EulerTerm::VX] = 2.0*_uInf - u;
     _dataGhostState[EulerTerm::VY] = 2.0*_vInf - v;
     _dataGhostState[EulerTerm::VZ] = 2.0*_wInf - w;
     _dataGhostState[EulerTerm::P] = 2.0*_pressure - pInnerState;
     _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
					  _dataGhostState[EulerTerm::RHO]);
     _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
					  _dataGhostState[EulerTerm::VX] +
					  _dataGhostState[EulerTerm::VY]*
					  _dataGhostState[EulerTerm::VY] +
					  _dataGhostState[EulerTerm::VZ]*
					  _dataGhostState[EulerTerm::VZ]);

     _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				      + 0.5*_dataGhostState[EulerTerm::RHO]*
				      _dataGhostState[EulerTerm::V]*
				      _dataGhostState[EulerTerm::V])/
       _dataGhostState[EulerTerm::RHO];
     
     _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*R);
     
     // set the ghost state starting from the physical data
     _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
  }

  if(fabs(machInner) < 1.0)
  {
    // subsonic outlet case
    if (RMinRmax > 0.)
    {

    CFreal vnInner = _dataInnerState[EulerTerm::VX]*nx + _dataInnerState[EulerTerm::VY] * ny + _dataInnerState[EulerTerm::VZ] * nz;
    CFreal Ubound = _dataInnerState[EulerTerm::VX] + (((RPlus_n + RMin_n)*0.5) - vnInner)*nx;
    CFreal Vbound = _dataInnerState[EulerTerm::VY] + (((RPlus_n + RMin_n)*0.5) - vnInner)*ny;
    CFreal Wbound = _dataInnerState[EulerTerm::VZ] + (((RPlus_n + RMin_n)*0.5) - vnInner)*nz;

     _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
     _dataGhostState[EulerTerm::VX] = 2.0*Ubound - u;
     _dataGhostState[EulerTerm::VY] = 2.0*Vbound - v;
     _dataGhostState[EulerTerm::VZ] = 2.0*Wbound - w;
     _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
					  _dataGhostState[EulerTerm::VX] +
					  _dataGhostState[EulerTerm::VY]*
					  _dataGhostState[EulerTerm::VY] +
					  _dataGhostState[EulerTerm::VZ]*
					  _dataGhostState[EulerTerm::VZ]);
     _dataGhostState[EulerTerm::P] = _dataInnerState[EulerTerm::P];
     _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				      + 0.5*_dataGhostState[EulerTerm::RHO]*
				      _dataGhostState[EulerTerm::V]*
				      _dataGhostState[EulerTerm::V])/
       _dataGhostState[EulerTerm::RHO];
     _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
					  _dataGhostState[EulerTerm::RHO]);
     
     _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*R);
     _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
    }

   // subsonic inlet case
    if (RMinRmax < 0.)
    {

    CFreal  Ubound = _uInf + (((RPlus_n + RMin_n)*0.5) - vnInf)*nx;
    CFreal  Vbound = _vInf + (((RPlus_n + RMin_n)*0.5) - vnInf)*ny;
    CFreal  Wbound = _wInf + (((RPlus_n + RMin_n)*0.5) - vnInf)*nz;

     _dataGhostState[EulerTerm::RHO] = _pressure/(R*_temperature);
     _dataGhostState[EulerTerm::VX] = 2.0*Ubound - u;
     _dataGhostState[EulerTerm::VY] = 2.0*Vbound - v;
     _dataGhostState[EulerTerm::VZ] = 2.0*Wbound - w;
     _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
					  _dataGhostState[EulerTerm::VX] +
					  _dataGhostState[EulerTerm::VY]*
					  _dataGhostState[EulerTerm::VY] +
					  _dataGhostState[EulerTerm::VZ]*
					  _dataGhostState[EulerTerm::VZ]);
     _dataGhostState[EulerTerm::P] = _pressure;
     _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				      + 0.5*_dataGhostState[EulerTerm::RHO]*
				      _dataGhostState[EulerTerm::V]*
				      _dataGhostState[EulerTerm::V])/
       _dataGhostState[EulerTerm::RHO];

     _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
					  _dataGhostState[EulerTerm::RHO]);
     
     _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*R);
     _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
    }
  }
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
