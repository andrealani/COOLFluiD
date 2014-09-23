#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DPhotosphere.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"

#include <gsl/gsl_sf_legendre.h>

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

MethodCommandProvider<MirrorMHD3DPhotosphere, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DPhotosphereFVMCCProvider("MirrorMHD3DPhotosphereFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DPhotosphere::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("rhoFixed", "Non-dimensional density value that is to be fixed in the ghost cells.");
  options.addConfigOption< CFreal >("TFixed", "Dimensional temperature value that is to be fixed in the ghost cells.");
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DPhotosphere::MirrorMHD3DPhotosphere(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  socket_Almreal("Almreal"),
  socket_Almimg("Almimg"),
  socket_Blmreal("Blmreal"),
  socket_Blmimg("Blmimg"),
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

  _TFixed = 1.5e6;
  setParameter("TFixed",&_TFixed);
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DPhotosphere::~MirrorMHD3DPhotosphere()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DPhotosphere::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DPhotosphere::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
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

//////////////////////////////////////////////////////////////////////////////

void MirrorMHD3DPhotosphere::computePFSSMagneticField(const RealVector& stateCoordsSpherical,
                                        RealVector& BPFSSSpherical)
{
  const CFreal r = stateCoordsSpherical[0];
  const CFreal theta = stateCoordsSpherical[1];
  const CFreal phi = stateCoordsSpherical[2];

  DataHandle<std::vector<CFreal> > Almreal  = socket_Almreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Almimg  = socket_Almimg.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmreal  = socket_Blmreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmimg  = socket_Blmimg.getDataHandle();

  const CFuint nbLModes = _varSet->getNbLModes();

  CFreal Br = 0.0, Btheta = 0.0, Bphi = 0.0;

  for (CFuint l = 0; l <= nbLModes; ++l) {
      for (CFuint m = 0; m <= l; ++m) {
          // Assign the real part of B(r,phi)          
          const CFreal Ylmampl = gsl_sf_legendre_sphPlm(l,m,cos(theta));
          const CFreal Brreal = -(Almreal[l])[m]*(double)l*pow(r,(double)l-1.0)+(Blmreal[l])[m]*((double)l+1.0)*pow(r,-(double)l-2.0);
          const CFreal Brimg = -(Almimg[l])[m]*(double)l*pow(r,(double)l-1.0)+(Blmimg[l])[m]*((double)l+1.0)*pow(r,-(double)l-2.0);
          const CFreal BthetaCoeff = -1.0/(r*sin(theta));
          const CFreal Bphireal = BthetaCoeff*(double)m*((Almreal[l])[m]*pow(r,(double)l)+(Blmreal[l])[m]*pow(r,-(double)l-1.0));
          const CFreal Bphiimg = BthetaCoeff*(double)m*((Almimg[l])[m]*pow(r,(double)l)+(Blmimg[l])[m]*pow(r,-(double)l-1.0));

          const CFreal Brampl = sqrt(Brreal*Brreal+Brimg*Brimg);
          const CFreal Bphiampl = sqrt(Bphireal*Bphireal+Bphiimg*Bphiimg);

          // TBD: In order to obtain the correct polarity of the magnetic dipoles on the photosphere -= is used instead of += for Br, Bphi and Btheta. This should be checked.
          const CFreal Brangle = atan2(Brimg,Brreal);
          Br -= Brampl*Ylmampl*cos(Brangle+(double)m*phi);

          const CFreal Bphiangle = atan2(Bphiimg,Bphireal);
          Bphi -= Bphiampl*Ylmampl*cos(Bphiangle+(double)m*phi+0.5*MathTools::MathConsts::CFrealPi());
      }
  }

  for (CFuint l = 1; l <= (nbLModes-1); ++l) {
      for (CFuint m = 0; m <= l; ++m) {
          const CFreal Rlm = sqrt(((double)l*(double)l-(double)m*(double)m)/(4.0*(double)l*(double)l-1.0));
          const CFreal Rlpl1m = sqrt((((double)l+1.0)*((double)l+1.0)-(double)m*(double)m)/(4.0*((double)l+1.0)*((double)l+1.0)-1.0));
          // Assign the real part of B(theta)          
          const CFreal Ylmampl = gsl_sf_legendre_sphPlm(l,m,cos(theta));
          const CFreal BthetaCoeff = -1.0/(r*sin(theta));
          const CFreal Bthetareal = BthetaCoeff*(Rlm*((double)l-1.0)*((Almreal[l-1])[m]*pow(r,(double)l-1.0)+(Blmreal[l-1])[m]*pow(r,-(double)l))
                                            -Rlpl1m*((double)l+2.0)*((Almreal[l+1])[m]*pow(r,(double)l+1.0)+(Blmreal[l+1])[m]*pow(r,-(double)l-2.0)));
          const CFreal Bthetaimg = BthetaCoeff*(Rlm*((double)l-1.0)*((Almimg[l-1])[m]*pow(r,(double)l-1.0)+(Blmimg[l-1])[m]*pow(r,-(double)l))
                                            -Rlpl1m*((double)l+2.0)*((Almimg[l+1])[m]*pow(r,(double)l+1.0)+(Blmimg[l+1])[m]*pow(r,-(double)l-2.0)));

          const CFreal Bthetaampl = sqrt(Bthetareal*Bthetareal+Bthetaimg*Bthetaimg);

          const CFreal Bthetaangle = atan2(Bthetaimg,Bthetareal);
          Btheta -= Bthetaampl*Ylmampl*cos(Bthetaangle+(double)m*phi);
      }
  }

  BPFSSSpherical[0] = Br;
  BPFSSSpherical[1] = Btheta;
  BPFSSSpherical[2] = Bphi;
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DPhotosphere::setGhostState(GeometricEntity *const face)
{
  DataHandle<std::vector<CFreal> > Almreal = socket_Almreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Almimg = socket_Almimg.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmreal = socket_Blmreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmimg = socket_Blmimg.getDataHandle();

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  RealVector innerStateCoordsSpherical(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector ghostStateCoordsSpherical(Framework::PhysicalModelStack::getActive()->getDim());

  //const CFreal lRef = _varSet->getLRef();

  const RealVector innerStateCoords = innerState->getCoordinates()/**lRef*/;
  const RealVector ghostStateCoords = ghostState->getCoordinates()/**lRef*/;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  // set the transformation matrices between Cartesian and spherical coordinate systems
  _varSet->setTransformationMatrices(innerStateCoords,innerStateCoordsSpherical,_cartesianSphericalTMInnerState,_sphericalCartesianTMInnerState);
  _varSet->setTransformationMatrices(ghostStateCoords,ghostStateCoordsSpherical,_cartesianSphericalTMGhostState,_sphericalCartesianTMGhostState);

  RealVector BCartesianInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector BSphericalInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector BCartesianGhostState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector BSphericalGhostState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VCartesianInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VSphericalInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VCartesianGhostState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VSphericalGhostState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector BPFSSGhostState(Framework::PhysicalModelStack::getActive()->getDim());

  BCartesianInnerState[0] = _dataInnerState[MHDTerm::BX];
  BCartesianInnerState[1] = _dataInnerState[MHDTerm::BY];
  BCartesianInnerState[2] = _dataInnerState[MHDTerm::BZ];

  VCartesianInnerState[0] = _dataInnerState[MHDTerm::VX];
  VCartesianInnerState[1] = _dataInnerState[MHDTerm::VY];
  VCartesianInnerState[2] = _dataInnerState[MHDTerm::VZ];

  BSphericalInnerState = _cartesianSphericalTMInnerState*BCartesianInnerState;
  VSphericalInnerState = _cartesianSphericalTMInnerState*VCartesianInnerState;

  VSphericalGhostState[0] = 0.0;
  VSphericalGhostState[1] = VSphericalInnerState[1];
  VSphericalGhostState[2] = -VSphericalInnerState[2];

  VCartesianGhostState = _sphericalCartesianTMGhostState*VSphericalGhostState;

  computePFSSMagneticField(ghostStateCoordsSpherical,BPFSSGhostState);

  BSphericalGhostState[0] = BPFSSGhostState[0]; 
  BSphericalGhostState[1] = BSphericalInnerState[1];
  BSphericalGhostState[2] = BSphericalInnerState[2];

  BCartesianGhostState = _sphericalCartesianTMGhostState*BSphericalGhostState;

  const CFreal VelCartesianGhostState = sqrt(VCartesianGhostState[0]*VCartesianGhostState[0]
			      + VCartesianGhostState[1]*VCartesianGhostState[1]
			      + VCartesianGhostState[2]*VCartesianGhostState[2]);

  const CFreal BMagCartesianGhostState = sqrt(BCartesianGhostState[0]*BCartesianGhostState[0]
                              + BCartesianGhostState[1]*BCartesianGhostState[1]
                              + BCartesianGhostState[2]*BCartesianGhostState[2]);

  // magnetic permeability at vacuum
  const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();

  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  // reference values used to non-dimensionalize the ideal MHD equations are assumed to occur at this boundary
  const CFreal nRef = _varSet->getNRef();
  const CFreal BRef = _varSet->getBRef();

  const CFreal rhoRef = nRef*(mp+me);
  const CFreal vRef = BRef / sqrt(mu0*rhoRef);
  
  const CFreal pFixed = (nRef*2.0*k*_TFixed)/(rhoRef*vRef*vRef); 

  // all boundary conditions are according to Carla Jacobs's suggestions
  _dataGhostState[MHDTerm::RHO] = 2.*_rhoFixed - _dataInnerState[MHDTerm::RHO];
  _dataGhostState[MHDTerm::VX] = VCartesianGhostState[0];
  _dataGhostState[MHDTerm::VY] = VCartesianGhostState[1];
  _dataGhostState[MHDTerm::VZ] = VCartesianGhostState[2];
  _dataGhostState[MHDTerm::BX] = BCartesianGhostState[0]/*/BRef*/; 
  _dataGhostState[MHDTerm::BY] = BCartesianGhostState[1]/*/BRef*/;
  _dataGhostState[MHDTerm::BZ] = BCartesianGhostState[2]/*/BRef*/;
  _dataGhostState[MHDTerm::V] = VelCartesianGhostState;
  _dataGhostState[MHDTerm::P] = 2.*pFixed - _dataInnerState[MHDTerm::P]; 
  _dataGhostState[MHDTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDTerm::P]/_dataGhostState[MHDTerm::RHO]);
  _dataGhostState[MHDTerm::B] = BMagCartesianGhostState/*/BRef*/;
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > MirrorMHD3DPhotosphere::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();
  result.push_back(&socket_Almreal);
  result.push_back(&socket_Almimg);
  result.push_back(&socket_Blmreal);
  result.push_back(&socket_Blmimg);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
