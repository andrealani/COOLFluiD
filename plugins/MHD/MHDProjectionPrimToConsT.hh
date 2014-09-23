#ifndef COOLFluiD_Physics_MHD_MHDProjectionPrimToConsT_hh
#define COOLFluiD_Physics_MHD_MHDProjectionPrimToConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformerT.hh"
#include "MHD/MHD2DProjectionPrimT.hh"
#include "MHD/MHD2DProjectionConsT.hh"
#include "MHD/MHD3DProjectionPrimT.hh"
#include "MHD/MHD3DProjectionConsT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive 
 * to conservative variables in 2D
 *
 * @author Andrea Lani
 */
template <>
class VarSetTransformerT<Physics::MHD::MHD2DProjectionPrimT, 
			 Physics::MHD::MHD2DProjectionConsT, 
			 NOTYPE> {
public:
  
  /// Constructor
  HOST_DEVICE VarSetTransformerT(Physics::MHD::MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
  
  /// Default constructor
  HOST_DEVICE VarSetTransformerT() {}
    
  /// Destructor
  HOST_DEVICE virtual ~VarSetTransformerT() {}
  
  /// set the model data
  HOST_DEVICE void setModelData(Physics::MHD::MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) {m_dco = dco;}
  
  /// Transform a state into another one
  HOST_DEVICE void transform(const CFreal *const state, CFreal *const result)
  {
    using namespace COOLFluiD::Physics::MHD;
    
    const CFreal gammaMinus1 = m_dco->gamma - 1.;
    const CFreal rho = state[0];
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal w = state[3];
    const CFreal Bx = state[4];
    const CFreal By = state[5];
    const CFreal Bz = state[6];
    const CFreal p = state[7];
    const CFreal phi = state[8];
    const CFreal rhoV2 = rho*(u*u + v*v + w*w);
    const CFreal B2 = (Bx*Bx + By*By + Bz*Bz);
    
    result[0] = rho;
    result[1] = rho*u;
    result[2] = rho*v;
    result[3] = rho*w;
    result[4] = Bx;
    result[5] = By;
    result[6] = Bz;
    result[7] = (p/gammaMinus1) + 0.5*(rhoV2+B2);
    result[8] = phi;
  }
  
  /// Transform a state into another one from reference precomputed
  HOST_DEVICE void transformFromRef(const CFreal *const data, CFreal *const result)
  {
    using namespace COOLFluiD::Physics::MHD;
    
    const CFreal gammaMinus1 = m_dco->gamma - 1.;
    
    const CFreal rhoV2 = data[MHDProjectionTerm::RHO]*
      (data[MHDProjectionTerm::VX]*data[MHDProjectionTerm::VX] + 
       data[MHDProjectionTerm::VY]*data[MHDProjectionTerm::VY] + 
       data[MHDProjectionTerm::VZ]*data[MHDProjectionTerm::VZ]);
    
    const CFreal B2 = data[MHDProjectionTerm::BX]*data[MHDProjectionTerm::BX] + 
      data[MHDProjectionTerm::BY]*data[MHDProjectionTerm::BY] + 
      data[MHDProjectionTerm::BZ]*data[MHDProjectionTerm::BZ];
    
    result[0] = data[MHDProjectionTerm::RHO];
    result[1] = data[MHDProjectionTerm::RHO]*data[MHDProjectionTerm::VX];
    result[2] = data[MHDProjectionTerm::RHO]*data[MHDProjectionTerm::VY];
    result[3] = data[MHDProjectionTerm::RHO]*data[MHDProjectionTerm::VZ];
    result[4] = data[MHDProjectionTerm::BX];
    result[5] = data[MHDProjectionTerm::BY];
    result[6] = data[MHDProjectionTerm::BZ];
    result[7] = (data[MHDProjectionTerm::P]/gammaMinus1) + 0.5*(rhoV2+B2);
    result[8] = data[MHDProjectionTerm::PHI];
  }
  
private:
  
  /// pointer to the model data
  Physics::MHD::MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* m_dco;
  
}; // end of class MHDProjectionPrimToConsT

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive 
 * to conservative variables in 3D
 *
 * @author Andrea Lani
 */
template <>
class VarSetTransformerT< Physics::MHD::MHD3DProjectionPrimT, Physics::MHD::MHD3DProjectionConsT, NOTYPE> :  
      public VarSetTransformerT< Physics::MHD::MHD2DProjectionPrimT, Physics::MHD::MHD2DProjectionConsT, NOTYPE> {
  
public:
  
  /// Constructor
  HOST_DEVICE VarSetTransformerT(Physics::MHD::MHDProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    VarSetTransformerT< Physics::MHD::MHD2DProjectionPrimT, Physics::MHD::MHD2DProjectionConsT, NOTYPE>(dco) {}
  
  /// Default constructor
  HOST_DEVICE VarSetTransformerT() : 
    VarSetTransformerT< Physics::MHD::MHD2DProjectionPrimT, Physics::MHD::MHD2DProjectionConsT, NOTYPE>() {}
  
  /// Destructor
  HOST_DEVICE ~VarSetTransformerT() {}
  
}; // end of class MHDProjectionPrimToConsT
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionPrimToConsT_hh
