#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DConsToRhoiViTiT_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DConsToRhoiViTiT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "MultiFluidMHD/EulerMFMHD2DRhoiViTiT.hh"
#include "MultiFluidMHD/EulerMFMHD2DConsT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a GPU-enabled transformer of variables from conservative variable to primitive 
 * [rhoi vi ui Ti]
 * AAL: This one is NOT used when you select RhoiViTi as vars
 *
 * @author Alejandro Alvarez
 * @author Isaac Alonso
 *
 */
template <>
class VarSetTransformerT<Physics::MultiFluidMHD::EulerMFMHD2DConsT, 
			 Physics::MultiFluidMHD::EulerMFMHD2DRhoiViTiT, 
			 NOTYPE> {

public:
  
  /// Constructor
  HOST_DEVICE VarSetTransformerT(Physics::MultiFluidMHD::EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
  
  /// Default constructor
  HOST_DEVICE VarSetTransformerT() {}
    
  /// Destructor
  HOST_DEVICE virtual ~VarSetTransformerT() {}
  
  /// set the model data
  HOST_DEVICE void setModelData(Physics::MultiFluidMHD::EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* dco) {m_dco = dco;}
  
    /// Transform a state into another one
  HOST_DEVICE void transform(const CFreal *const state, CFreal *const result)
  {
    const CFuint nbSpecies = m_dco->nbSpecies;
    const CFuint nbMomentum = m_dco->nbMomentum;
    const CFuint nbEnergyEqs = m_dco->nbEnergyEqs;
    const CFuint endEM = 8;
  
    //Electro Magnetic Field Needs no tranformation   
    result[0] = state[0];
    result[1] = state[1];  
    result[2] = state[2];
    result[3] = state[3];
    result[4] = state[4];
    result[5] = state[5]; 
    result[6] = state[6];
    result[7] = state[7];
  
  
    //The densities continue wthout transformation
  
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      result[endEM + ie] = state[endEM + ie];
    }
  
    //Set the velocities Ui = RhoiUi/Rhoi and Vi = RhoiVi/Rhoi
  
    for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
      result[endEM + nbSpecies + 2*ie] = state[endEM + nbSpecies + 2*ie]/state[endEM + ie];
      result[endEM + nbSpecies + 2*ie + 1] = state[endEM + nbSpecies + 2*ie +1]/state[endEM + ie];
    }
  
    const bool isLeake = m_dco->isLeake;

    if (isLeake){
      //ions
      //set the energy parameters
      const CFreal gamma = m_dco->gamma;	// gamma = 5/3
      const CFreal K_B = m_dco->K;		// Boltzmann constant
    
      const CFreal m_p = m_dco->molecularMass3;
    
      const CFreal R_i = K_B/m_p;				// ions gas constant
      const CFreal R_p = 2*R_i;				// Plasma gas constant (ions + electrons)
      const CFreal Cp_p = (gamma/(gamma-1))*R_p;	
      const CFreal Cv_p = (1/(gamma-1))*R_p;
      const CFreal u_i = result[endEM + nbSpecies];
      const CFreal v_i = result[endEM + nbSpecies + 1];
      const CFreal V2_i = u_i*u_i + v_i*v_i;
    
      result[endEM + nbSpecies + 2*nbSpecies] = (state[endEM + 3*nbSpecies]/state[endEM] - 0.5*V2_i)/Cv_p; 
    
      //neutrals
      //set the energy parameters
      const CFreal m_n = m_dco->molecularMass2;
    
      const CFreal R_n = K_B/m_n;				// neutrals gas constant
      const CFreal Cp_n = (gamma/(gamma-1))*R_n;	
      const CFreal Cv_n = (1/(gamma-1))*R_n;
      const CFreal u_n = result[endEM + nbSpecies + 2];
      const CFreal v_n = result[endEM + nbSpecies + 3];
      const CFreal V2_n = u_n*u_n + v_n*v_n;
    
      result[endEM + 3*nbSpecies + 1] = (state[endEM + 3*nbSpecies + 1]/state[endEM + 1] - 0.5*V2_n)/Cv_n;    
    
    }
  
    else{   
      //set the energy parameters
      const CFreal gamma = m_dco->gamma;
      const CFreal K_gas = m_dco->K;
      const CFreal m_e = m_dco->molecularMass1;
      const CFreal m_n = m_dco->molecularMass2;
      const CFreal m_p = m_dco->molecularMass3; 
    
      //set the molar masses of the species (should be changed in the future)

      _m_i[0] = m_e;
      _m_i[1] = m_n;
      _m_i[2] = m_p;  
    
      //Set Ti = (RhoiEi/Rhoi - 0.5*V²)/Cv
    
      for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    
        const CFreal u = state[endEM + nbSpecies + 2*ie]/state[endEM + ie];
        const CFreal v = state[endEM + nbSpecies + 2*ie +1]/state[endEM + ie];
        const CFreal V2 = u*u + v*v;
        const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
        const CFreal R_gas = K_gas/_m_i[ie];
        const CFreal c_v = c_p - R_gas;  
      
        result[endEM + 3*nbSpecies + ie] = (state[endEM + 3*nbSpecies + ie]/state[endEM + ie] - 0.5*V2)/c_v; 
      }
    }
  }



  /// Transform a state into another one from reference precomputed
  HOST_DEVICE void transformFromRef(const CFreal *const data, CFreal *const result)
  {
    const CFuint nbSpecies = m_dco->nbSpecies;
    const CFuint nbMomentum   = m_dco->nbMomentum;
    const CFuint nbEnergyEqs  = m_dco->nbEnergyEqs;
    const CFuint endEM = 8;
 
    //Electro Magnetic Field Needs no tranformation   
    result[0] = data[0];
    result[1] = data[1];  
    result[2] = data[2];
    result[3] = data[3];
    result[4] = data[4];
    result[5] = data[5]; 
    result[6] = data[6];
    result[7] = data[7];  
  
    //Density of the mixture
    const CFreal rho = data[8];   //PTERM::RHO
    //The densities Rhoi = Rho*yi
  
    const CFuint firstSpecies = m_dco->firstSpecies;
  
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      result[endEM + ie] = data[firstSpecies + ie]*rho;
    }
  
    //Set the momentum Ui and Vi
  
    const CFuint firstVelocity = m_dco->firstVelocity;  
    for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
      result[endEM + nbSpecies + 2*ie] = data[firstVelocity + 2*ie];
      result[endEM + nbSpecies + 2*ie + 1] = data[firstVelocity + 2*ie +1];
    }
    
  
    //Set Ti 
    const CFuint firstTemperature = m_dco->firstTemperature;  
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
  
      result[endEM + 3*nbSpecies + ie] = data[firstTemperature + 4*ie]; 
    }
  }
  
private:

  Physics::MultiFluidMHD::EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* m_dco;
  CFreal _m_i[3];

}; // end of class EulerMFMHD2DConsToRhoiViTiT

//////////////////////////////////////////////////////////////////////////////


  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DConsToRhoiViTiT_hh
