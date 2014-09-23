#ifndef COOLFluiD_Physics_MHD_MHD2DProjectionPrimT_hh
#define COOLFluiD_Physics_MHD_MHD2DProjectionPrimT_hh

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
class MHD2DProjectionPrimT : public MHD2DProjectionVarSetT {
public: // function
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD2DProjectionPrimT(MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    MHD2DProjectionVarSetT(dco) {}
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  HOST_DEVICE MHD2DProjectionPrimT() : MHD2DProjectionVarSetT() {}
  
  /**
   * Default destructor
   */
  HOST_DEVICE ~MHD2DProjectionPrimT() {}
  
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
   const CFreal V2 = u*u + v*v + w*w;
   const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
   const CFreal gamma = this->m_dco->gamma;
   
   data[MHDProjectionTerm::VX] = u;
   data[MHDProjectionTerm::VY] = v;
   data[MHDProjectionTerm::VZ] = w;
   data[MHDProjectionTerm::BX] = Bx;
   data[MHDProjectionTerm::BY] = By;
   data[MHDProjectionTerm::BZ] = Bz;
   data[MHDProjectionTerm::RHO] = rho;
   data[MHDProjectionTerm::V] = sqrt(V2);
   data[MHDProjectionTerm::B] = sqrt(B2);
   data[MHDProjectionTerm::P] = state[7];
   data[MHDProjectionTerm::PHI] = state[8];
   data[MHDProjectionTerm::A] = sqrt(gamma*data[MHDProjectionTerm::P]/rho);
   data[MHDProjectionTerm::GAMMA] = gamma;
   data[MHDTerm::XP] =  node[XX];
   data[MHDTerm::YP] =  node[YY];
  }
  
}; // end of class MHD2DProjectionPrimT

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DProjectionPrimT_hh
