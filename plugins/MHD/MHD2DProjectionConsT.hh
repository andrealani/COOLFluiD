#ifndef COOLFluiD_Physics_MHD_MHD2DProjectionConsT_hh
#define COOLFluiD_Physics_MHD_MHD2DProjectionConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD2DProjectionVarSetT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 2D for projection scheme
 * for conservative variables
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD2DProjectionConsT : public MHD2DProjectionVarSetT {
public: // function
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD2DProjectionConsT(MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    MHD2DProjectionVarSetT(dco) {}
 
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD2DProjectionConsT() : MHD2DProjectionVarSetT() {}
  
  /**
   * Default destructor
   */
  HOST_DEVICE ~MHD2DProjectionConsT() {}
  
  /// Compute the physical data starting ffrom the corresponding state variables
  HOST_DEVICE void computePhysicalData(CFreal* state, CFreal* node, CFreal* data) 
  {
    const CFreal rho = state[0];
    const CFreal u = state[1]/rho;
    const CFreal v = state[2]/rho;
    const CFreal w = state[3]/rho;
    const CFreal Bx = state[4];
    const CFreal By = state[5];
    const CFreal Bz = state[6];
    const CFreal phi = state[8];
    const CFreal V2 = u*u + v*v;
    const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
    const CFreal gamma = this->m_dco->gamma;
    
    data[MHDTerm::VX] = u;
    data[MHDTerm::VY] = v;
    data[MHDTerm::VZ] = w;
    data[MHDTerm::BX] = Bx;
    data[MHDTerm::BY] = By;
    data[MHDTerm::BZ] = Bz;
    data[MHDTerm::B] = sqrt(B2);
    data[MHDTerm::RHO] = rho;
    data[MHDTerm::V] = sqrt(V2);
    data[MHDTerm::P] = (gamma - 1.)*(state[7] - 0.5*(rho*V2 + B2));
    data[MHDProjectionTerm::PHI] = phi;
    data[MHDTerm::A] = sqrt(gamma*data[MHDTerm::P]/rho);
    data[MHDProjectionTerm::GAMMA] = gamma;
    data[MHDTerm::XP] =  node[XX];
    data[MHDTerm::YP] =  node[YY];
  }
  
}; // end of class MHD2DProjectionConsT

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DProjectionConsT_hh
