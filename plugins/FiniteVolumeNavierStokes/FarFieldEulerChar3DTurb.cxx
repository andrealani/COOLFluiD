#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/FarFieldEulerChar3DTurb.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/MeshData.hh"
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

MethodCommandProvider<FarFieldEulerChar3DTurb, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
farFieldEulerChar3DTurbFVMCCProvider("FarFieldEulerChar3DTurbFVMCC");

//////////////////////////////////////////////////////////////////////////////

void FarFieldEulerChar3DTurb::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("TurbVars","Freestream K, Omega values (assigned by defaulty to litterature values)");
}

/////////////////////////////////////////////////////////////////////////////

FarFieldEulerChar3DTurb::FarFieldEulerChar3DTurb(const std::string& name) :
  FarFieldEulerChar3D(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL)
{

   addConfigOptionsTo(this);

   _turbVars = vector<CFreal>();
   setParameter("TurbVars",&_turbVars);

}

//////////////////////////////////////////////////////////////////////////////

FarFieldEulerChar3DTurb::~FarFieldEulerChar3DTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void FarFieldEulerChar3DTurb::setup()
{
  FarFieldEulerChar3D::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();
  cf_assert(_varSetTurb.isNotNull());

  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb3DVarSet>();
  cf_assert(_diffVarTurb.isNotNull());

  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  if(_turbVars.size() == 0){
    _turbVars.resize(_varSetTurb->getModel()->getNbScalarVars(0));
  }

  //This is for one equation model
  if(_turbVars.size() == 1){
    throw Common::NotImplementedException
      (FromHere(), "FarFieldEulerChar3DTurb not implemented for SA model");
  }

  // for k-Omega model
  else{
    cf_assert(_turbVars.size() == 2);

    ///@todo this is only valid for k-Omega model
    //check the adimensionalisation - velocities are in (m/s), temperature in (Kelvin)

    //Values taken from: F.R. Menter: Two-Eq Eddy-Viscosity Turbulence Models for Engineering Applications (Aug 1994)
    const CFreal L = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
    const CFreal Uinf = sqrt(_uInf*_uInf + _vInf*_vInf);
    //Wilcox BC
    _turbVars[1] = Uinf/L;
    //Menter BC
    //_turbVars[1] = 10.*Uinf/L;

    const CFreal pdim = _pressure * _varSetTurb->getModel()->getPressRef();
    const CFreal Tdim = _temperature * _varSetTurb->getModel()->getTempRef();

    const CFreal muInf = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
                            (_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
    
    //upper bound
    const CFreal muTurbInf = muInf/100.;
    //lower bound
    //const CFreal muTurbInf = muInf/100000.;

    _turbVars[0] = _turbVars[1] * muTurbInf ;
  }

  CFout << "Setting the turbulent values to :     k = " << _turbVars[0] << 
         "\n                                  Omega = " << _turbVars[1] << "\n";

  //Check that the initial values for the turbulent variables have been set
  cf_assert(_turbVars.size() == _varSetTurb->getModel()->getNbScalarVars(0));

  const CFuint firstScalarVar = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbScalarVars = _varSetTurb->getModel()->getNbScalarVars(0);

  const RealVector& refValues = _varSetTurb->getModel()->getReferencePhysicalData();
  for(CFuint iVar=0; iVar < nbScalarVars ;iVar++)
  {
    _turbVars[iVar] /= refValues[firstScalarVar + iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void FarFieldEulerChar3DTurb::setGhostState(GeometricEntity *const face)
{
   State *const innerState = face->getState(0);
   State *const ghostState = face->getState(1);

   const CFreal gamma = _varSetTurb->getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.0;
   const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
   const CFreal R = _varSetTurb->getModel()->getR();

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
   _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

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

  /// supersonic outlet ---------------------------------------------------------------------------
  if (machInner >= 1.0){
    const CFuint nbVars = (*innerState).size(); // returns 7 for 3D-turb computations
    for(CFuint iVar = 0; iVar < nbVars; iVar++){
      (*ghostState)[iVar] = (*innerState)[iVar];
    }
  }

  /// supersonic inlet ----------------------------------------------------------------------------
  if (machInner <= -1.0){
    const CFreal rhoInf = _pressure/(R*_temperature);

    // set all the physical data corresponding to the ghost state
    _dataGhostState[EulerTerm::RHO] = 2.0*rhoInf - rho;
    _dataGhostState[EulerTerm::VX] = 2.0*_uInf - u;
    _dataGhostState[EulerTerm::VY] = 2.0*_vInf - v;
    _dataGhostState[EulerTerm::VZ] = 2.0*_wInf - w;
    _dataGhostState[EulerTerm::V] = sqrt( _dataGhostState[EulerTerm::VX]*_dataGhostState[EulerTerm::VX]
                                        + _dataGhostState[EulerTerm::VY]*_dataGhostState[EulerTerm::VY]
                                        + _dataGhostState[EulerTerm::VZ]*_dataGhostState[EulerTerm::VZ]);
    _dataGhostState[EulerTerm::P] = 2.0*_pressure - pInnerState;
    _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);
    _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				      + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*_dataGhostState[EulerTerm::V])/
                    _dataGhostState[EulerTerm::RHO];
    
    _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*R);
    
    const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
    const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
    for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++){
      _dataGhostState[iK + iTurb] = 2.0*_turbVars[iTurb] - _dataInnerState[iK + iTurb];
    }

    // set the ghost state starting from the physical data
    _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
  }

  if(fabs(machInner) < 1.0){
    /// subsonic outlet ---------------------------------------------------------------------------
    if (RMinRmax > 0.){
      CFreal vnInner = _dataInnerState[EulerTerm::VX]*nx + _dataInnerState[EulerTerm::VY] * ny + _dataInnerState[EulerTerm::VZ] * nz;
      CFreal Ubound = _dataInnerState[EulerTerm::VX] + (((RPlus_n + RMin_n)*0.5) - vnInner)*nx;
      CFreal Vbound = _dataInnerState[EulerTerm::VY] + (((RPlus_n + RMin_n)*0.5) - vnInner)*ny;
      CFreal Wbound = _dataInnerState[EulerTerm::VZ] + (((RPlus_n + RMin_n)*0.5) - vnInner)*nz;

      _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
      _dataGhostState[EulerTerm::VX] = 2.0*Ubound - u;
      _dataGhostState[EulerTerm::VY] = 2.0*Vbound - v;
      _dataGhostState[EulerTerm::VZ] = 2.0*Wbound - w;
      _dataGhostState[EulerTerm::V] = sqrt( _dataGhostState[EulerTerm::VX]*_dataGhostState[EulerTerm::VX]
                                          + _dataGhostState[EulerTerm::VY]*_dataGhostState[EulerTerm::VY]
                                          + _dataGhostState[EulerTerm::VZ]*_dataGhostState[EulerTerm::VZ]);
      _dataGhostState[EulerTerm::P] = _dataInnerState[EulerTerm::P];
      _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
	 			      + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*_dataGhostState[EulerTerm::V])/
                    _dataGhostState[EulerTerm::RHO];
      _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);
      _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*R);
      
      const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
      const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
      for(CFuint iTurb = iK; iTurb < iK + nbTurbVars; iTurb++){
        _dataGhostState[iTurb] = _dataInnerState[iTurb];
      }
      _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
    }

    /// subsonic inlet ----------------------------------------------------------------------------
    if (RMinRmax < 0.){
      CFreal  Ubound = _uInf + (((RPlus_n + RMin_n)*0.5) - vnInf)*nx;
      CFreal  Vbound = _vInf + (((RPlus_n + RMin_n)*0.5) - vnInf)*ny;
      CFreal  Wbound = _wInf + (((RPlus_n + RMin_n)*0.5) - vnInf)*nz;

      _dataGhostState[EulerTerm::RHO] = _pressure/(R*_temperature);
      _dataGhostState[EulerTerm::VX] = 2.0*Ubound - u;
      _dataGhostState[EulerTerm::VY] = 2.0*Vbound - v;
      _dataGhostState[EulerTerm::VZ] = 2.0*Wbound - w;
      _dataGhostState[EulerTerm::V] = sqrt( _dataGhostState[EulerTerm::VX]*_dataGhostState[EulerTerm::VX]
                                          + _dataGhostState[EulerTerm::VY]*_dataGhostState[EulerTerm::VY]
                                          + _dataGhostState[EulerTerm::VZ]*_dataGhostState[EulerTerm::VZ]);
      _dataGhostState[EulerTerm::P] = _pressure;
      _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
	 	  		      + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*_dataGhostState[EulerTerm::V])/
                      _dataGhostState[EulerTerm::RHO];

      _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);
      _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*R);
      
      const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
      const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
      for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++){
        _dataGhostState[iK + iTurb] = 2.0*_turbVars[iTurb] - _dataInnerState[iK + iTurb];
      }

      _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
