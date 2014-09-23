#include "MHD3DProjectionPolytropicVarSet.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MeshData.hh"

#ifdef CF_HAVE_GSL
#include <gsl/gsl_sf_legendre.h>
#endif

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicVarSet::MHD3DProjectionPolytropicVarSet(Common::SafePtr<Framework::BaseTerm> term) :
  ConvectiveVarSet(term),
  _model(term.d_castTo<MHDProjectionPolytropicTerm>()),
  _BDipole(),
  _tanakaFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicVarSet::~MHD3DProjectionPolytropicVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicVarSet::computeFlux (const RealVector& data,
					 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal rho = data[MHDProjectionPolytropicTerm::RHO];
  const CFreal u = data[MHDProjectionPolytropicTerm::VX];
  const CFreal v = data[MHDProjectionPolytropicTerm::VY];
  const CFreal w = data[MHDProjectionPolytropicTerm::VZ];
  const CFreal Bx = data[MHDProjectionPolytropicTerm::BX];
  const CFreal By = data[MHDProjectionPolytropicTerm::BY];
  const CFreal Bz = data[MHDProjectionPolytropicTerm::BZ];
  const CFreal Vn = u*nx  + v*ny + w*nz;
  const CFreal Bn = Bx*nx + By*ny + Bz*nz;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal p = data[MHDProjectionPolytropicTerm::P];
  const CFreal phi = data[MHDProjectionPolytropicTerm::PHI];
  const CFreal P = p + 0.5*B2;

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _fluxArray[0] = Vn*rho;
  _fluxArray[1] = Vn*rho*u - Bn*Bx + P*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + P*ny;
  _fluxArray[3] = Vn*rho*w - Bn*Bz + P*nz;
  _fluxArray[4] = (v*Bx - By*u)*ny + (w*Bx - Bz*u)*nz + phi*nx;
  _fluxArray[5] = (u*By - Bx*v)*nx + (w*By - Bz*v)*nz + phi*ny;
  _fluxArray[6] = (u*Bz - Bx*w)*nx + (v*Bz - By*w)*ny + phi*nz;
  _fluxArray[7] = refSpeedSq*Bn;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicVarSet::computeStateFlux (const RealVector& data)
{
  const CFreal rho = data[MHDProjectionPolytropicTerm::RHO];
  const CFreal u = data[MHDProjectionPolytropicTerm::VX];
  const CFreal v = data[MHDProjectionPolytropicTerm::VY];
  const CFreal w = data[MHDProjectionPolytropicTerm::VZ];
  const CFreal Bx = data[MHDProjectionPolytropicTerm::BX];
  const CFreal By = data[MHDProjectionPolytropicTerm::BY];
  const CFreal Bz = data[MHDProjectionPolytropicTerm::BZ];
  const CFreal BxBy = Bx*By;
  const CFreal BxBz = Bx*Bz;
  const CFreal ByBz = By*Bz;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal p = data[MHDProjectionPolytropicTerm::P];
  const CFreal phi = data[MHDProjectionPolytropicTerm::PHI];
  const CFreal P = p + 0.5*B2;
  const CFreal rhoU = rho*u;
  const CFreal rhoV = rho*v;
  const CFreal rhoW = rho*w;
  const CFreal rhoUU = rhoU*u;
  const CFreal rhoVV = rhoV*v;
  const CFreal rhoWW = rhoW*w;
  const CFreal rhoUV = rhoU*v;
  const CFreal rhoUW = rhoU*w;
  const CFreal rhoVW = rhoV*w;
  const CFreal vBx = v*Bx;
  const CFreal wBx = w*Bx;
  const CFreal uBy = u*By;
  const CFreal wBy = w*By;
  const CFreal uBz = u*Bz;
  const CFreal vBz = v*Bz;

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _physFlux(0,XX) = rhoU;
  _physFlux(0,YY) = rhoV;
  _physFlux(0,ZZ) = rhoW;

  _physFlux(1,XX) = rhoUU + Bx*Bx + P;
  _physFlux(1,YY) = rhoUV + BxBy;
  _physFlux(1,ZZ) = rhoUW + BxBz;

  _physFlux(2,XX) = rhoUV + BxBy;
  _physFlux(2,YY) = rhoVV + By*By+ P;
  _physFlux(2,ZZ) = rhoVW + ByBz;

  _physFlux(3,XX) = rhoUW + BxBz;
  _physFlux(3,YY) = rhoVW + ByBz;
  _physFlux(3,ZZ) = rhoWW + Bz*Bz + P;

  _physFlux(4,XX) = phi;
  _physFlux(4,YY) = vBx - uBy;
  _physFlux(4,ZZ) = wBx - uBz;

  _physFlux(5,XX) = uBy - vBx;
  _physFlux(5,YY) = phi;
  _physFlux(5,ZZ) = wBy - vBz;

  _physFlux(6,XX) = uBz - wBx;
  _physFlux(6,YY) = vBz - wBy;
  _physFlux(6,ZZ) = phi;

  _physFlux(7,XX) = refSpeedSq*Bx;
  _physFlux(7,YY) = refSpeedSq*By;
  _physFlux(7,ZZ) = refSpeedSq*Bz;
}

//////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicVarSet::setTransformationMatrices(const RealVector& coords,
                           RealVector& coordsSpherical,
			   RealMatrix& carSphTransMat,
                           RealMatrix& sphCarTransMat)
{
        const CFreal x = coords[0];
        const CFreal y = coords[1];
        const CFreal z = coords[2];
        const CFreal r = sqrt(x*x+y*y+z*z);
        const CFreal R = sqrt(x*x+y*y);

	coordsSpherical[0] = r;
	coordsSpherical[1] = atan2(R,z);
        coordsSpherical[2] = atan2(y,x);

        carSphTransMat(0,0) = x/r;
        carSphTransMat(1,0) = (x*z)/(R*r);
        carSphTransMat(2,0) = -y/R;

        carSphTransMat(0,1) = y/r;
        carSphTransMat(1,1) = (y*z)/(R*r);
        carSphTransMat(2,1) = x/R;

        carSphTransMat(0,2) = z/r;
        carSphTransMat(1,2) = -R/r;
        carSphTransMat(2,2) = 0.0;

        sphCarTransMat(0,0) = x/r;
        sphCarTransMat(1,0) = y/r;
        sphCarTransMat(2,0) = z/r;

        sphCarTransMat(0,1) = (x*z)/(R*r);
        sphCarTransMat(1,1) = (y*z)/(R*r);
        sphCarTransMat(2,1) = -R/r;

        sphCarTransMat(0,2) = -y/R;
        sphCarTransMat(1,2) = x/R;
        sphCarTransMat(2,2) = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD3DProjectionPolytropicVarSet::computeTanakaFluxPowell99Formulation(const RealVector& data,
                                                                        const RealVector& normals)
{
  // According to the explanations in Powell, K.G. et. al., JCP, Vol.154, pp.284-309, 1999

  const std::string potentialBType = getModel()->getPotentialBType();

  CFreal Bx0 = 0., By0 = 0., Bz0 = 0.;
  if (potentialBType == "Dipole") {
     computeMagneticDipole(data[MHDProjectionPolytropicTerm::XP],data[MHDProjectionPolytropicTerm::YP],data[MHDProjectionPolytropicTerm::ZP]);
     Bx0 = _BDipole[0];
     By0 = _BDipole[1];
     Bz0 = _BDipole[2];
  }
  if (potentialBType == "PFSS") {
     std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
     std::string datahandleName2 = nsp + "_BPFSS";

     DataHandle<std::vector<CFreal> > BPFSS = MeshDataStack::getActive()->getDataStorage()->getData<std::vector<CFreal> >(datahandleName2);

     // Boltzmann constant
     const CFreal k = 1.3806503e-23;

     // magnetic permeability at vacuum
     const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();

     // mass of proton
     const CFreal mp = 1.67262158e-27;

     // mass of electron
     const CFreal me = 9.10938188e-31;

     const CFreal nRef = getModel()->getNRef();
     const CFreal TRef = getModel()->getTRef();

     const CFreal rhoRef = nRef*(mp+me);
     const CFreal vRef = sqrt(2.0*k*TRef/mp);

     const CFreal BRef = sqrt(mu0*rhoRef*vRef*vRef)*1.e4;

     Bx0 = (BPFSS[_entityID])[0]/BRef;
     By0 = (BPFSS[_entityID])[1]/BRef;
     Bz0 = (BPFSS[_entityID])[2]/BRef;
  }

  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal rho = data[MHDProjectionPolytropicTerm::RHO];
  const CFreal u = data[MHDProjectionPolytropicTerm::VX];
  const CFreal v = data[MHDProjectionPolytropicTerm::VY];
  const CFreal w = data[MHDProjectionPolytropicTerm::VZ];
  const CFreal Bx1 = data[MHDProjectionPolytropicTerm::BX];
  const CFreal By1 = data[MHDProjectionPolytropicTerm::BY];
  const CFreal Bz1 = data[MHDProjectionPolytropicTerm::BZ];

  // n stands for normal to the face of the control volume

  const CFreal Vn = u*nx + v*ny + w*nz;
  const CFreal Bn0 = Bx0*nx + By0*ny + Bz0*nz;

  const CFreal Bn1 = Bx1*nx + By1*ny + Bz1*nz;
  const CFreal sqB1 = Bx1*Bx1 + By1*By1 + Bz1*Bz1;
  const CFreal B1dotB0 = Bx1*Bx0 + By1*By0 + Bz1*Bz0;
  const CFreal p = data[MHDProjectionPolytropicTerm::P];

  const CFreal P = p + 0.5*sqB1 + B1dotB0;

  const CFreal phi = data[MHDProjectionTerm::PHI];

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _tanakaFlux[0] = Vn*rho;
  _tanakaFlux[1] = Vn*rho*u - Bn1*Bx1 - Bn0*Bx1 - Bn1*Bx0 + P*nx;
  _tanakaFlux[2] = Vn*rho*v - Bn1*By1 - Bn0*By1 - Bn1*By0 + P*ny;
  _tanakaFlux[3] = Vn*rho*w - Bn1*Bz1 - Bn0*Bz1 - Bn1*Bz0 + P*nz;
  _tanakaFlux[4] = (v*Bx1 - By1*u)*ny + (w*Bx1 - Bz1*u)*nz + (v*Bx0 - By0*u)*ny + (w*Bx0 - Bz0*u)*nz + phi*nx;
  _tanakaFlux[5] = (u*By1 - Bx1*v)*nx + (w*By1 - Bz1*v)*nz + (u*By0 - Bx0*v)*nx + (w*By0 - Bz0*v)*nz + phi*ny;
  _tanakaFlux[6] = (u*Bz1 - Bx1*w)*nx + (v*Bz1 - By1*w)*ny + (u*Bz0 - Bx0*w)*nx + (v*Bz0 - By0*w)*ny + phi*nz;
  _tanakaFlux[7] = refSpeedSq*Bn1;

  return _tanakaFlux;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicVarSet::computeMagneticDipole(CFreal xCoord, CFreal yCoord, CFreal zCoord)
{
  // this formula is valid for a dipole at the origin
  const CFreal mX = getMX();
  const CFreal mY = getMY();
  const CFreal mZ = getMZ();

  const CFreal r = sqrt(xCoord*xCoord + yCoord*yCoord + zCoord*zCoord);
  const CFreal rCubeInv = 1.0/(r*r*r);
  const CFreal rSqInv = 1.0/(r*r);
  const CFreal mDotr = mX*xCoord + mY*yCoord + mZ*zCoord;

  const CFreal BxDipole = rCubeInv*(3.0*xCoord*rSqInv*mDotr - mX);
  const CFreal ByDipole = rCubeInv*(3.0*yCoord*rSqInv*mDotr - mY);
  const CFreal BzDipole = rCubeInv*(3.0*zCoord*rSqInv*mDotr - mZ);

  _BDipole[0] = BxDipole;
  _BDipole[1] = ByDipole;
  _BDipole[2] = BzDipole;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicVarSet::computePFSSMagneticField(const RealVector& stateCoordsSpherical,
                                        RealVector& BPFSSCartesian,
                                        RealMatrix& sphCarTransMat)
{
#ifdef CF_HAVE_GSL
  CFreal r = stateCoordsSpherical[0];
  const CFreal theta = stateCoordsSpherical[1];
  const CFreal phi = stateCoordsSpherical[2];
  const CFreal rTemp = r;

  const CFreal rSource = getModel()->getRSource();

  if (r >= rSource) {
     r = rSource;
  }

  std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName1, datahandleName2, datahandleName3, datahandleName4;

  const bool interpolationFlag = getModel()->getInterpolationFlag();

  if (interpolationFlag == true) {
	datahandleName1 = nsp + "_AlmrealEnd";
        datahandleName2 = nsp + "_AlmimgEnd";
        datahandleName3 = nsp + "_BlmrealEnd";
        datahandleName4 = nsp + "_BlmimgEnd";
  }

  if (interpolationFlag == false) {
  	datahandleName1 = nsp + "_AlmrealBegin";
  	datahandleName2 = nsp + "_AlmimgBegin";
  	datahandleName3 = nsp + "_BlmrealBegin";
  	datahandleName4 = nsp + "_BlmimgBegin";
  }

  DataHandle<std::vector<CFreal> > Almreal = MeshDataStack::getActive()->getDataStorage()->getData<std::vector<CFreal> >(datahandleName1);
  DataHandle<std::vector<CFreal> > Almimg = MeshDataStack::getActive()->getDataStorage()->getData<std::vector<CFreal> >(datahandleName2);
  DataHandle<std::vector<CFreal> > Blmreal = MeshDataStack::getActive()->getDataStorage()->getData<std::vector<CFreal> >(datahandleName3);
  DataHandle<std::vector<CFreal> > Blmimg = MeshDataStack::getActive()->getDataStorage()->getData<std::vector<CFreal> >(datahandleName4);

  const CFuint nbLModes = getModel()->getNbLModes();

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

  RealVector BPFSSSpherical(PhysicalModelStack::getActive()->getDim());

  BPFSSSpherical[0] = Br;
  if (rTemp < rSource) {
        BPFSSSpherical[1] = Btheta;
        BPFSSSpherical[2] = Bphi;
  } else {
        if (rTemp > rSource) {
              // B field is radial and decreases by r^2 beyond the source surface
              BPFSSSpherical[0] = Br/(rTemp*rTemp);
        }
        BPFSSSpherical[1] = 0.0;
        BPFSSSpherical[2] = 0.0;
  }

  BPFSSCartesian = sphCarTransMat*BPFSSSpherical;
#endif
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionPolytropicVarSet::getMaxEigenValue(const RealVector& data,
					       const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDProjectionPolytropicTerm::RHO];
  const CFreal B = data[MHDProjectionPolytropicTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDProjectionPolytropicTerm::P];
  const CFreal Vn = data[MHDProjectionPolytropicTerm::VX]*normal[XX] + data[MHDProjectionPolytropicTerm::VY]*normal[YY] +
    data[MHDProjectionPolytropicTerm::VZ]*normal[ZZ];
  const CFreal Bn = data[MHDProjectionPolytropicTerm::BX]*normal[XX] + data[MHDProjectionPolytropicTerm::BY]*normal[YY] +
    data[MHDProjectionPolytropicTerm::BZ]*normal[ZZ];
  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
				  - 4.0*getModel()->getGamma()*p*Bn*Bn*
				  invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal maxEigenValue = max(refSpeed,(Vn + cf));
  return maxEigenValue;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionPolytropicVarSet::getMaxAbsEigenValue(const RealVector& data,
						  const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDProjectionPolytropicTerm::RHO];
  const CFreal B = data[MHDProjectionPolytropicTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDProjectionPolytropicTerm::P];
  const CFreal Vn = data[MHDProjectionPolytropicTerm::VX]*normal[XX] + data[MHDProjectionPolytropicTerm::VY]*normal[YY] +
    data[MHDProjectionPolytropicTerm::VZ]*normal[ZZ];
  const CFreal Bn = data[MHDProjectionPolytropicTerm::BX]*normal[XX] + data[MHDProjectionPolytropicTerm::BY]*normal[YY] +
    data[MHDProjectionPolytropicTerm::BZ]*normal[ZZ];
  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
          - 4.0*getModel()->getGamma()*p*Bn*Bn*
          invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal maxEigenValue = max(refSpeed,(abs(Vn) + cf));
  return maxEigenValue;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicVarSet::computeEigenValues(const RealVector& data,
                                               const RealVector& normal,
                                               RealVector& result)
{
  // eigenvalues for the flows not involving any background potential magnetic field are also correctly
  // calculated in this function since B0 field will automatically be zero in that case
  // eigenvalues based on total magnetic field, Btotal

  const std::string potentialBType = getModel()->getPotentialBType();

  CFreal Bx0 = 0., By0 = 0., Bz0 = 0.;
  if (potentialBType == "Dipole") {
     computeMagneticDipole(data[MHDProjectionPolytropicTerm::XP],data[MHDProjectionPolytropicTerm::YP],data[MHDProjectionPolytropicTerm::ZP]);
     Bx0 = _BDipole[0];
     By0 = _BDipole[1];
     Bz0 = _BDipole[2];
  }
  if (potentialBType == "PFSS") {
     std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
     std::string datahandleName2 = nsp + "_BPFSS";

     DataHandle<std::vector<CFreal> > BPFSS = MeshDataStack::getActive()->getDataStorage()->getData<std::vector<CFreal> >(datahandleName2);          

     // Boltzmann constant
     const CFreal k = 1.3806503e-23;
     
     // magnetic permeability at vacuum
     const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();
  
     // mass of proton
     const CFreal mp = 1.67262158e-27;
  
     // mass of electron
     const CFreal me = 9.10938188e-31;
  
     const CFreal nRef = getModel()->getNRef();
     const CFreal TRef = getModel()->getTRef();
  
     const CFreal rhoRef = nRef*(mp+me);
     const CFreal vRef = sqrt(2.0*k*TRef/mp);

     const CFreal BRef = sqrt(mu0*rhoRef*vRef*vRef)*1.e4;

     Bx0 = (BPFSS[_entityID])[0]/BRef;
     By0 = (BPFSS[_entityID])[1]/BRef;
     Bz0 = (BPFSS[_entityID])[2]/BRef;
  }

  const CFreal Bx1 = data[MHDProjectionPolytropicTerm::BX];
  const CFreal By1 = data[MHDProjectionPolytropicTerm::BY];
  const CFreal Bz1 = data[MHDProjectionPolytropicTerm::BZ];
  const CFreal BxTotal = Bx1 + Bx0;
  const CFreal ByTotal = By1 + By0;
  const CFreal BzTotal = Bz1 + Bz0;

  const CFreal invRho = 1./data[MHDProjectionPolytropicTerm::RHO];
  //const CFreal B = data[MHDProjectionPolytropicTerm::B];
  const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
  //const CFreal B2 = B*B;
  const CFreal p = data[MHDProjectionPolytropicTerm::P];
  const CFreal Vn = data[MHDProjectionPolytropicTerm::VX]*normal[XX] + data[MHDProjectionPolytropicTerm::VY]*normal[YY] +
    data[MHDProjectionPolytropicTerm::VZ]*normal[ZZ];
  //const CFreal Bn = data[MHDProjectionPolytropicTerm::BX]*normal[XX] + data[MHDProjectionPolytropicTerm::BY]*normal[YY] +
  //data[MHDProjectionPolytropicTerm::BZ]*normal[ZZ];
  const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
  const CFreal sqrbar = sqrt(data[MHDProjectionPolytropicTerm::RHO]);
  const CFreal astar2 = (getModel()->getGamma()*p + sqBTotal)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
                                  - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho));
  CFreal cs2 = 0.5*(astar2 - sqrt(astar2*astar2
                                  - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho));

  const CFreal ca = std::fabs(BnTotal)/sqrbar;
  const CFreal cf = sqrt(cf2);
  const CFreal cs = sqrt(cs2);

  const CFreal refSpeed = getModel()->getRefSpeed();

  result[0] = Vn - cf;
  result[1] = Vn - ca;
  result[2] = Vn - cs;
  //result[3] = Vn;
  result[3] = refSpeed;
  result[4] = -refSpeed;
  result[5] = Vn + cs;
  result[6] = Vn + ca;
  result[7] = Vn + cf;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
