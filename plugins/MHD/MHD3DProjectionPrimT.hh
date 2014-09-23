#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionPrimT_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionPrimT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD3DProjectionVarSetT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for projection scheme
 * for conservative variables
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD3DProjectionPrimT : public MHD3DProjectionVarSetT {
public: // function
 
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD3DProjectionPrimT(MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    MHD3DProjectionVarSetT(dco) {}
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD3DProjectionPrimT() : MHD3DProjectionVarSetT() {}
  
  /**
   * Default destructor
   */
  HOST_DEVICE ~MHD3DProjectionPrimT() {}
  
  /// Compute the physical data starting ffrom the corresponding state variables
  HOST_DEVICE void computePhysicalData(CFreal* state, CFreal* node, CFreal* data) 
  {
    const CFreal rho = state[0];
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal w = state[3];
    const CFreal Bx = state[4];
    const CFreal By = state[5];
    const CFreal Bz = state[6];
    const CFreal phi = state[8];
    const CFreal V2 = u*u + v*v + w*w;
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
    data[MHDTerm::P] = state[7];
    data[MHDProjectionTerm::PHI] = phi;
    data[MHDTerm::A] = sqrt(gamma*data[MHDTerm::P]/rho);
    data[MHDTerm::GAMMA] = gamma;
    data[MHDTerm::XP] =  node[XX];
    data[MHDTerm::YP] =  node[YY];
    data[MHDTerm::ZP] =  node[ZZ];
  }
  
}; // end of class MHD3DProjectionPrimT

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionPrimT_hh
