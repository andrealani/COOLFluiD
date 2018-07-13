#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DHalfRhoiViTiT_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DHalfRhoiViTiT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MultiFluidMHD2DHalfTwoSpeciesVarSetT.hh"
#include "Maxwell/Maxwell2DProjectionVarSetT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a GPU-enabled multi-fluid Euler physical 
 * model 2.5D for Rhoi, Vi, Ti
 * variables:
 * Bx, By, Bz, Ex, Ey, Ez, Psi, Phi,
 * {rhoi}, {Ui, Vi, Wi}, {Ti}
 * 
 * PHYSICAL DATA ARRAY:
 * Bx, By, Bz, Ex, Ey, Ez, Psi, Phi,
 * rho, XP, YP, ZP,
 * {yi},
 * {Ui, Vi, Wi},
 * {Ti, pi, ai, Hi} 
 * 
 * 
 * @author Alejandro Alvarez
 * @author Isaac Alonso
 * 
 */
class EulerMFMHD2DHalfRhoiViTiT : public MultiFluidMHD2DHalfTwoSpeciesVarSetT<Maxwell::Maxwell2DProjectionVarSetT> {
public: // classes

  typedef EulerMFMHDTerm PTERM;
  
  /**
   * Constructor
   * @see MultiFluidMHDModel
   */
  
   /**
   * Constructor
   */
  HOST_DEVICE EulerMFMHD2DHalfRhoiViTiT(EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    MultiFluidMHD2DHalfTwoSpeciesVarSetT(dco) {}
 
  /**
   * Constructor
   */
  HOST_DEVICE EulerMFMHD2DHalfRhoiViTiT() : MultiFluidMHD2DHalfTwoSpeciesVarSetT() {}


  /**
   * Default destructor
   */
  HOST_DEVICE ~EulerMFMHD2DHalfRhoiViTiT() {};
  
  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  //virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  //CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  //virtual void computeJacobians();

  /**
   * Set the jacobian matrix
   */
  //virtual void computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob);
  
  /**
   * Split the jacobian
   */
  //void splitJacobian(RealMatrix& jacobPlus,
  //		     RealMatrix& jacobMin,
  //		     RealVector& eValues,
  //		     const RealVector& normal);
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  //void computeEigenValuesVectors(RealMatrix& rightEv,
  //				 RealMatrix& leftEv,
  //				 RealVector& eValues,
  //				 const RealVector& normal);
  
  /**
   * Get the speed
   */
  //CFreal getSpeed(const Framework::State& state) const;
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  //void setDimensionalValues(const Framework::State& state, RealVector& result);
  
  /**
   * Give adimensional values to the dimensional state variables
   */
  //void setAdimensionalValues(const Framework::State& state, RealVector& result);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  HOST_DEVICE void computePhysicalData(CFreal* state, CFreal* data)
  {
    //printf("EulerMFMHD2DHalfRhoiViTiT::computePhysicalData \n"); 
    const CFuint nbSpecies    = m_dco->nbSpecies;
    const CFuint endEM = 8;
    const CFuint firstSpecies = m_dco->firstSpecies;  
    const CFuint firstVelocity = m_dco->firstVelocity;   
    const CFuint firstTemperature = m_dco->firstTemperature;

    const CFreal m_e = m_dco->molecularMass1;
    const CFreal m_n = m_dco->molecularMass2;
    const CFreal m_p = m_dco->molecularMass3; 


    const CFreal gamma = m_dco->gamma;
    const CFreal K_gas = m_dco->K;


    //DEBUG
    //printf("gamma %e \t K_gas %e \t m_e %e \t m_n %e \t m_p %e \n",
    //        gamma,      K_gas,      m_e,      m_n,      m_p);


    //printf("nbSpecies %d \t firstSpecies %d \t firstVelocity %d \t firstTemperature %d \n",
    //        nbSpecies,      firstSpecies,      firstVelocity,      firstTemperature);


    //set the molar masses of the species (should be changed in the future)
    _m_i[0] = m_e;
    _m_i[1] = m_n;
    _m_i[2] = m_p;
    
    data[PTERM::BX] = state[0];
    data[PTERM::BY] = state[1]; 
    data[PTERM::BZ] = state[2];  
    data[PTERM::EX] = state[3];  
    data[PTERM::EY] = state[4];  
    data[PTERM::EZ] = state[5];  
    data[PTERM::PSI] = state[6];  
    data[PTERM::PHI] = state[7];  
    
    const bool isLeake = m_dco->isLeake;

    // plasma + neutrals model
    if(isLeake){
    
      //Total density
      CFreal rho = 0.0;
      CFreal rho_i = state[endEM];
      CFreal rho_n = state[endEM + 1];
    
      rho = (1 + m_e/m_p)*rho_i + rho_n;  // rho_Total = rho_i + rho_e + rho_n = (1 + m_e/m_i)rho_i + rho_n
      data[PTERM::RHO] = rho;
      const CFreal ovRho = 1./rho;
    
      cf_assert(data[PTERM::RHO] > 0.);
    
      //Partial densities
      data[firstSpecies] =  rho_i*ovRho;
      data[firstSpecies + 1] =  rho_n*ovRho;  


      //Velocities
      const CFreal u_i = state[endEM + nbSpecies];	//Ions x-velocity
      const CFreal v_i = state[endEM + nbSpecies + 1];	//Ions y-velocity
      const CFreal w_i = state[endEM + nbSpecies + 2];	//Ions z-velocity
      data[firstVelocity]     = u_i;
      data[firstVelocity + 1] = v_i; 
      data[firstVelocity + 2] = w_i;
    
      const CFreal u_n = state[endEM + nbSpecies + 3];	//Neutrals x-velocity
      const CFreal v_n = state[endEM + nbSpecies + 4];	//Neutrals y-velocity
      const CFreal w_n = state[endEM + nbSpecies + 5];	//Neutrals z-velocity
      data[firstVelocity + 3] = u_n;
      data[firstVelocity + 4] = v_n; 
      data[firstVelocity + 5] = w_n;
    
      //Energy Variables: Ti, pi, ai, Hi
      
      const CFreal gamma = m_dco->gamma;	// gamma = 5/3
      const CFreal K_B = m_dco->K;     // Boltzmann constant
              

      //plasma
      const CFuint dim = 3;
      const CFreal V2_i = u_i*u_i + v_i*v_i + w_i*w_i;
      const CFreal T_i = state[endEM + nbSpecies + dim*nbSpecies];       
      
      const CFreal R_i = K_B/m_p;				// ions gas constant
      const CFreal R_p = 2.*R_i;				// Plasma gas constant (ions + electrons)
      const CFreal Cv_p = R_p/(gamma-1.);
      const CFreal Cp_p = gamma*Cv_p;
        
      data[firstTemperature] = T_i;			//Temperature
      data[firstTemperature + 1] = T_i*R_p*rho_i;		//pressure (pe + pi)
      data[firstTemperature + 2] = sqrt(gamma*R_p*T_i);	//sound speed (pe + pi)
      data[firstTemperature + 3] = 0.5*V2_i + Cp_p*T_i;	//total enthaply of species (pe + pi) 
      
      //neutrals
      const CFreal V2_n = u_n*u_n + v_n*v_n + w_n*w_n;
      const CFreal T_n = state[endEM + nbSpecies + dim*nbSpecies + 1];       
  
      const CFreal R_n  = K_B/m_n;             // neutrals gas constant
      const CFreal Cv_n = R_n/(gamma-1.);      // Cv for neutrals 
      const CFreal Cp_n = gamma*Cv_n;          // Cp for neutrals
        
      data[firstTemperature + 4] = T_n;				//Temperature
      data[firstTemperature + 4 + 1] = T_n*R_n*rho_n;		//pressure
      data[firstTemperature + 4 + 2] = sqrt(gamma*R_n*T_n);	//sound speed
      data[firstTemperature + 4 + 3] = 0.5*V2_n + Cp_n*T_n;	//total enthalpy of species i     
    }
  
    else{
 
      //set the total density
      CFreal rho = 0.0;
      for (CFuint ie = 0; ie < nbSpecies; ++ie) {
        rho += state[endEM + ie];
      }

      data[PTERM::RHO] = rho;
      const CFreal ovRho = 1./rho;
    
      //set the energy parameters

      
      //set the molar masses of the species (should be changed in the future)
      _m_i[0] = m_e;
      _m_i[1] = m_n;
      _m_i[2] = m_p;
    
      //printf("nbSpecies %d \n", nbSpecies);
      //set the species mass fraction
      for (CFuint ie = 0; ie < nbSpecies; ++ie) {
        const CFreal rhoi = state[endEM + ie];
        
        data[firstSpecies + ie] = rhoi*ovRho;
        //printf("data[%d] %e \n", firstSpecies + ie, data[firstSpecies + ie]);
        //set the species velocities in 2.5D 
        const CFuint dim = 3;
        const CFreal ui = state[endEM + nbSpecies + dim*ie];
        const CFreal vi = state[endEM + nbSpecies + dim*ie + 1];
        const CFreal wi = state[endEM + nbSpecies + dim*ie + 2];
        data[firstVelocity + dim*ie] = ui;
        data[firstVelocity + dim*ie + 1] = vi;
        data[firstVelocity + dim*ie + 2] = wi;
      
        
      
      
        const CFreal mi = _m_i[ie];    
        const CFreal V2 = ui*ui + vi*vi + wi*wi;
        const CFreal Ti = state[endEM + nbSpecies + dim*nbSpecies + ie];       
    
        const CFreal c_p = (gamma/(gamma-1))*(K_gas/mi);
        const CFreal R_gas = K_gas/mi;
        //const CFreal c_v = c_p - R_gas;
        //printf("FirstTemperature %d \n", firstTemperature);
        data[firstTemperature + 4*ie] = Ti;//Temperature
        data[firstTemperature + 4*ie + 1] = Ti*R_gas*rhoi;//pressure
        data[firstTemperature + 4*ie + 2] = sqrt(gamma*R_gas*Ti);//sound speed
        data[firstTemperature + 4*ie + 3] = 0.5*V2 + c_p*Ti;//total enthaply of species i   
        //cout << "ie = "<< ie <<"\n";
        //cout << "V2 = "<< V2 <<"\n";
        //cout << "Ti = "<< Ti <<"\n";
      }
    }
  }  


  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  HOST_DEVICE void computeStateFromPhysicalData(const RealVector& data, Framework::State& state)
  {
    //printf("EulerMFMHD2DHalfRhoiViTi::computeStateFromPhysicalData \n"); 
    const CFuint nbSpecies = m_dco->nbSpecies;
    const CFuint endEM = 8;  
  
    state[0] = data[PTERM::BX];
    state[1] = data[PTERM::BY]; 
    state[2] = data[PTERM::BZ];  
    state[3] = data[PTERM::EX];  
    state[4] = data[PTERM::EY];  
    state[5] = data[PTERM::EZ];  
    state[6] = data[PTERM::PSI];  
    state[7] = data[PTERM::PHI];  
  
    const CFreal rho = data[PTERM::RHO];
  
    //set the species densities
    const CFuint firstSpecies = m_dco->firstSpecies;
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      state[endEM + ie] = data[firstSpecies + ie]*rho; //rhoi = Rho*yi
    }  
  
    //set the species velocities in 2.5D 
    const CFuint firstVelocity = m_dco->firstVelocity;
    const CFuint dim = 3;
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      state[endEM + nbSpecies + dim*ie] = data[firstVelocity + dim*ie];
      state[endEM + nbSpecies + dim*ie + 1] = data[firstVelocity + dim*ie + 1];
      state[endEM + nbSpecies + dim*ie + 2] = data[firstVelocity + dim*ie + 2];
    
    }

 
   //set the temperatures
    const CFuint firstTemperature = m_dco->firstTemperature;
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      state[endEM + nbSpecies + dim*nbSpecies + ie] = data[firstTemperature + 4*ie];
    }  
    //for (CFuint ie = 0; ie < endEM + 5*nbSpecies; ++ie) {
    //   printf("state[%d] = %f \n", ie, state[ie]);
    //} 
  }






  /// Compute the perturbed physical data
  //virtual void computePerturbedPhysicalData(const Framework::State& state,
  //					    const RealVector& pdataBkp,
  //					    RealVector& pdata,
  //					    CFuint iVar);
  
  /**
   * Returns true if the state fed doesn't have unphysical values
   * Overrides the standard definition in ConvectiveVarSet.
   */
  //bool isValid(const RealVector& data);

  

  
}; // end of class EulerMFMHD2DHalfRhoiViTi
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHD2DHalfRhoiViTi_hh
