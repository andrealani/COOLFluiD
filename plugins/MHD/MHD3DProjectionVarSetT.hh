#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionVarSetT_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionVarSetT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for projection scheme
 * for 3D MHD physical model.
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 */
class MHD3DProjectionVarSetT {
  
public: // classes
  
  enum {DIM=3, NBEQS=9, DATASIZE = 16};
  typedef MHDProjectionTerm PTERM;
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD3DProjectionVarSetT(MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) :
    m_dco(dco) {}
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD3DProjectionVarSetT() {}
  
  /**
   * Set the model data
   */
  HOST_DEVICE void setModelData(MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) {m_dco = dco;}
  
  /**
   * Default destructor
   */
  HOST_DEVICE virtual ~MHD3DProjectionVarSetT() {}
  
  /**
   * Get the x-component of the magnetic dipole moment (mX)
   */
  HOST_DEVICE  CFreal getMX() {return m_dco->mX;}
  
  /**
   * Get the y-component of the magnetic dipole moment (mY)
   */
  HOST_DEVICE  CFreal getMY() {return m_dco->mY;}
  
  /**
   * Get the z-component of the magnetic dipole moment (mZ)
   */
  HOST_DEVICE  CFreal getMZ() {return m_dco->mZ;}
  
  /// Computes the convective flux projected on a normal
  HOST_DEVICE  void getFlux(CFreal* data, CFreal* normal, CFreal* flux) 
  {
    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];
    const CFreal nz = normal[ZZ];
    const CFreal rho = data[MHDTerm::RHO];
    const CFreal u = data[MHDTerm::VX];
    const CFreal v = data[MHDTerm::VY];
    const CFreal w = data[MHDTerm::VZ];
    const CFreal Bx = data[MHDTerm::BX];
    const CFreal By = data[MHDTerm::BY];
    const CFreal Bz = data[MHDTerm::BZ];
    const CFreal Vn = u*nx  + v*ny + w*nz;
    const CFreal Bn = Bx*nx + By*ny + Bz*nz;
    const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
    const CFreal VdotB = u*Bx  + v*By + w*Bz;
    const CFreal p = data[MHDTerm::P];
    const CFreal phi = data[MHDProjectionTerm::PHI];
    const CFreal P = p + 0.5*B2;
    const CFreal gamma = m_dco->gamma;
    const CFreal E = p/(gamma - 1.0) + 0.5*(rho*(u*u + v*v + w*w) + B2);
    
    const CFreal refSpeed = m_dco->refSpeed;
    const CFreal refSpeedSq = refSpeed*refSpeed;
    
    flux[0] = Vn*rho;
    flux[1] = Vn*rho*u - Bn*Bx + P*nx;
    flux[2] = Vn*rho*v - Bn*By + P*ny;
    flux[3] = Vn*rho*w - Bn*Bz + P*nz;
    flux[4] = (v*Bx - By*u)*ny + (w*Bx - Bz*u)*nz + phi*nx;
    flux[5] = (u*By - Bx*v)*nx + (w*By - Bz*v)*nz + phi*ny;
    flux[6] = (u*Bz - Bx*w)*nx + (v*Bz - By*w)*ny + phi*nz;
    flux[7] = Vn*(E + P) - Bn*VdotB;
    flux[8] = refSpeedSq*Bn;
  }
    
  /// Compute the convective flux according to Powell, K.G. et. al. (JCP Vol.154 pp.284-309 1999)
  HOST_DEVICE void computeTanakaFluxPowell99Formulation(CFreal*  data, 
							CFreal* normals, 
							CFreal* flux)
  {
    CFreal Bx0 = 0.0, By0 = 0.0, Bz0 = 0.0;
    if (m_dco->potentialBType == MHDTerm::DIPOLE) {
      CFreal BDipole[3];
      computeMagneticDipole(data[MHDTerm::XP],data[MHDTerm::YP],data[MHDTerm::ZP], &BDipole[0]);
      Bx0 = BDipole[XX];
      By0 = BDipole[YY];
      Bz0 = BDipole[ZZ];
    }
    
    const CFreal nx = normals[XX];
    const CFreal ny = normals[YY];
    const CFreal nz = normals[ZZ];
    const CFreal rho = data[MHDTerm::RHO];
    const CFreal u = data[MHDTerm::VX];
    const CFreal v = data[MHDTerm::VY];
    const CFreal w = data[MHDTerm::VZ];
    const CFreal Bx1 = data[MHDTerm::BX];
    const CFreal By1 = data[MHDTerm::BY];
    const CFreal Bz1 = data[MHDTerm::BZ];
    
    // n stands for normal to the face of the control volume
    const CFreal Vn = u*nx + v*ny + w*nz;
    const CFreal Bn0 = Bx0*nx + By0*ny + Bz0*nz;
    const CFreal Bn1 = Bx1*nx + By1*ny + Bz1*nz;
    const CFreal sqB1 = Bx1*Bx1 + By1*By1 + Bz1*Bz1;
    const CFreal VdotB1 = u*Bx1  + v*By1 + w*Bz1;
    const CFreal V2 = u*u + v*v + w*w;
    const CFreal B1dotB0 = Bx1*Bx0 + By1*By0 + Bz1*Bz0;
    const CFreal p = data[MHDTerm::P];
    const CFreal P = p + 0.5*sqB1 + B1dotB0;
    const CFreal E1 = p/(m_dco->gamma - 1.) + 0.5 * (rho*V2 + sqB1);
    const CFreal phi = data[MHDProjectionTerm::PHI];
    const CFreal refSpeed = m_dco->refSpeed;
    const CFreal refSpeedSq = refSpeed*refSpeed;
    const CFreal H = E1 + p + 0.5*sqB1;
    
    flux[0] = Vn*rho;
    flux[1] = Vn*rho*u - Bn1*Bx1 - Bn0*Bx1 - Bn1*Bx0 + P*nx;
    flux[2] = Vn*rho*v - Bn1*By1 - Bn0*By1 - Bn1*By0 + P*ny;
    flux[3] = Vn*rho*w - Bn1*Bz1 - Bn0*Bz1 - Bn1*Bz0 + P*nz;
    flux[4] = (v*Bx1 - By1*u)*ny + (w*Bx1 - Bz1*u)*nz + (v*Bx0 - By0*u)*ny + (w*Bx0 - Bz0*u)*nz + phi*nx;
    flux[5] = (u*By1 - Bx1*v)*nx + (w*By1 - Bz1*v)*nz + (u*By0 - Bx0*v)*nx + (w*By0 - Bz0*v)*nz + phi*ny;
    flux[6] = (u*Bz1 - Bx1*w)*nx + (v*Bz1 - By1*w)*ny + (u*Bz0 - Bx0*w)*nx + (v*Bz0 - By0*w)*ny + phi*nz;
    flux[7] = Vn*H - Bn1*VdotB1 + Vn*B1dotB0 - Bn0*VdotB1;
    flux[8] = refSpeedSq*Bn1;
  }
  
  /// Set the vector of the eigenValues
  HOST_DEVICE  void computeEigenValues (CFreal* data, CFreal* normal, CFreal* eValues)
  {    
    // eigenvalues for the flows not involving any dipole field are also correctly
    // calculated in this function since B0 field will automatically be zero in that case
    // eigenvalues based on total magnetic field, Btotal
    CFreal Bx0 = 0.0, By0 = 0.0, Bz0 = 0.0;
    if (m_dco->potentialBType == MHDTerm::DIPOLE) {
      CFreal BDipole[3];
      computeMagneticDipole(data[MHDTerm::XP],data[MHDTerm::YP],data[MHDTerm::ZP], &BDipole[0]);
      Bx0 = BDipole[XX];
      By0 = BDipole[YY];
      Bz0 = BDipole[ZZ];
    }
    
    const CFreal refSpeed = m_dco->refSpeed;
    const CFreal gamma = m_dco->gamma;
    const CFreal BxTotal = data[MHDTerm::BX] + Bx0;
    const CFreal ByTotal = data[MHDTerm::BY] + By0;
    const CFreal BzTotal = data[MHDTerm::BZ] + Bz0;
    const CFreal invRho = 1./data[MHDTerm::RHO];
    const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
    const CFreal p = data[MHDTerm::P];
    const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];
    const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
    const CFreal sqrbar = sqrt(data[MHDTerm::RHO]);
    const CFreal astar2 = (gamma*p + sqBTotal)*invRho;
    const CFreal astarb = sqrt(astar2*astar2 - 4.0*gamma*p*BnTotal*BnTotal*invRho*invRho);
    const CFreal cf2 = 0.5*(astar2 + astarb);
    const CFreal cs2 = 0.5*(astar2 - astarb);
    const CFreal ca = fabs(BnTotal)/sqrbar;
    const CFreal cf = sqrt(cf2);
    const CFreal cs = sqrt(cs2);
    
    eValues[0] = Vn - cf;
    eValues[1] = Vn - ca;
    eValues[2] = Vn - cs;
    eValues[3] = Vn;
    eValues[4] = refSpeed;
    eValues[5] = -refSpeed;
    eValues[6] = Vn + cs;
    eValues[7] = Vn + ca;
    eValues[8] = Vn + cf;
  }

  /// Compute the magnetic dipole
  HOST_DEVICE  void computeMagneticDipole(CFreal xCoord, CFreal yCoord, 
					  CFreal zCoord, CFreal* BDipole)
  {
    // this formula is valid for a dipole at the origin
    const CFreal mX = getMX();
    const CFreal mY = getMY();
    const CFreal mZ = getMZ();
    
    const CFreal r = sqrt(xCoord*xCoord + yCoord*yCoord + zCoord*zCoord);
    const CFreal rCubeInv = 1.0/(r*r*r);
    const CFreal rSqInv = 1.0/(r*r);
    const CFreal mDotr = mX*xCoord + mY*yCoord + mZ*zCoord;
    const CFreal coeff = 3.0*rSqInv*mDotr;
    
    BDipole[XX] = rCubeInv*(coeff*xCoord - mX);
    BDipole[YY] = rCubeInv*(coeff*yCoord - mY);
    BDipole[ZZ] = rCubeInv*(coeff*zCoord - mZ);
  }
  
  /// Compute the maximum eigenvalue
  HOST_DEVICE  CFreal getMaxEigenValue(CFreal* data, CFreal* normal)
  {
    // eigenvalues for the flows not involving any dipole field are also correctly
    // calculated in this function since B0 field will automatically be zero in that case
    // eigenvalues based on total magnetic field, Btotal
    CFreal Bx0 = 0.0, By0 = 0.0, Bz0 = 0.0;
    if (m_dco->potentialBType == MHDTerm::DIPOLE) {
      CFreal BDipole[3];
      computeMagneticDipole(data[MHDTerm::XP],data[MHDTerm::YP],data[MHDTerm::ZP], &BDipole[0]);
      Bx0 = BDipole[XX];
      By0 = BDipole[YY];
      Bz0 = BDipole[ZZ];
    }
    
    const CFreal refSpeed = m_dco->refSpeed;
    const CFreal gamma = m_dco->gamma;
    const CFreal BxTotal = data[MHDTerm::BX] + Bx0;
    const CFreal ByTotal = data[MHDTerm::BY] + By0;
    const CFreal BzTotal = data[MHDTerm::BZ] + Bz0;
    const CFreal invRho = 1./data[MHDTerm::RHO];
    const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
    const CFreal p = data[MHDTerm::P];
    const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];
    const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
    const CFreal sqrbar = sqrt(data[MHDTerm::RHO]);
    const CFreal astar2 = (gamma*p + sqBTotal)*invRho;
    const CFreal astarb = sqrt(astar2*astar2 - 4.0*gamma*p*BnTotal*BnTotal*invRho*invRho);
    const CFreal cf2 = 0.5*(astar2 + astarb);
    const CFreal cf = sqrt(cf2);
    // const CFreal cf = sqrt(abs(cf2));
    const CFreal maxEigenValue = max(refSpeed,(Vn + cf));
    return maxEigenValue;
  }
  
protected:
  
  /// configurable options
  MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* m_dco;
  
}; // end of class MHD3DProjectionVarSetT

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionVarSetT_hh
