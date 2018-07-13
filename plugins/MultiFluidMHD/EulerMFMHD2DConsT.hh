#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DConsT_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MultiFluidMHD2DTwoSpeciesVarSetT.hh"
#include "Maxwell/Maxwell2DProjectionVarSetT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a multi-fluid (GPU-enabled)
 * Euler physical model 2D for conservative
 * variables:
 * Bx, By, Bz, Ex, Ey, Ez, Psi, Phi,
 * {rhoi},{rhoiUi, rhoiVi, RhoiW}, {rhoiEi} 
 * 
 * PHYSICAL DATA ARRAY:
 * Bx, By, Bz, Ex, Ey, Ez, Psi, Phi,
 * rho, XP, YP, ZP,
 * {yi},
 * {Ui, Vi, Wi},
 * {Ti, pi, ai, Hi} 
 * 
 * @author Alejandro Alvarez
 * @author Isaac Alonso
 */

class EulerMFMHD2DConsT : public MultiFluidMHD2DTwoSpeciesVarSetT<Maxwell::Maxwell2DProjectionVarSetT> { 
public:
  
  /**
   * Constructor
   * @see MultiFluidMHDModel
   */
  

   /**
   * Constructor
   */
  HOST_DEVICE EulerMFMHD2DConsT(EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* dco) : MultiFluidMHD2DTwoSpeciesVarSetT(dco) {};
 
  /**
   * Constructor
   */
  HOST_DEVICE EulerMFMHD2DConsT() : MultiFluidMHD2DTwoSpeciesVarSetT() {}; 


  /**
   * Default destructor
   */
  HOST_DEVICE ~EulerMFMHD2DConsT() {};
  
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
    //printf("EulerMFMHD2DConsT::computePhysicalData \n");

    const CFuint nbSpecies = m_dco->nbSpecies;
    //const CFuint nbMomentum = m_dco->nbMomentum;
    //const CFuint nbEnergyEqs  = m_dco->nbEnergyEqs;
    const CFuint endEM = 8;
    const CFuint firstSpecies = m_dco->firstSpecies;
    const CFuint firstVelocity = m_dco->firstVelocity;
    //const CFuint firstTemperature = m_dco->firstTemperature;
  
    data[PTERM::BX] = state[0];
    data[PTERM::BY] = state[1]; 
    data[PTERM::BZ] = state[2];  
    data[PTERM::EX] = state[3];  
    data[PTERM::EY] = state[4];  
    data[PTERM::EZ] = state[5];  
    data[PTERM::PSI] = state[6];  
    data[PTERM::PHI] = state[7];  
  
  
    //set the total density
    CFreal rho = 0.0;
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      rho += state[endEM + ie];
    }

    data[PTERM::RHO] = rho;
    const CFreal ovRho = 1./rho;
  
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
  
  
  //set the species mass fraction
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal rhoi = state[endEM + ie];
    data[firstSpecies + ie] = rhoi*ovRho;
    //set the species velocities in 2D
    const CFreal ui = state[endEM + nbSpecies + 2*ie]/rhoi;
    const CFreal vi = state[endEM + nbSpecies + 2*ie + 1]/rhoi;
    data[firstVelocity + 2*ie] = ui;
    data[firstVelocity + 2*ie + 1] = vi;
 
    //set the energy physical data
    const CFuint firstTemperature = m_dco->firstTemperature;
    
    const CFreal mi = _m_i[ie];    
    const CFreal V2 = ui*ui + vi*vi;
    const CFreal c_p = (gamma/(gamma-1))*(K_gas/mi);
    const CFreal R_gas = K_gas/mi;
    const CFreal c_v = c_p - R_gas;
    const CFreal Ti = (state[endEM + nbSpecies + 2*nbSpecies + ie] - rhoi*V2)/(rhoi*c_v);
    
    data[firstTemperature + 4*ie] = Ti;//Temperature
    data[firstTemperature + 4*ie + 1] = Ti*R_gas*rhoi;//pressure
    data[firstTemperature + 4*ie + 2] = sqrt(gamma*R_gas*Ti);//sound speed
    data[firstTemperature + 4*ie + 3] = 0.5*V2 + c_p*Ti;//total enthaply of species i   
  }    
}



  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  HOST_DEVICE void computeStateFromPhysicalData(CFreal* data, CFreal* state)
  {
    //printf("EulerMFMHD2DCons::computeStateFromPhysicalData \n");
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
  
    //set the species velocities in 2D 
    const CFuint firstVelocity = m_dco->firstVelocity;  
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      state[endEM + nbSpecies + 2*ie] = data[firstVelocity + 2*ie]*data[firstSpecies + ie]*rho;
      state[endEM + nbSpecies + 2*ie + 1] = data[firstVelocity + 2*ie + 1]*data[firstSpecies + ie]*rho;    
    }

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
 
    //set the Energies
    const CFuint firstTemperature = m_dco->firstTemperature;
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      const CFreal V2 = data[firstVelocity + 2*ie]*data[firstVelocity + 2*ie] + data[firstVelocity + 2*ie + 1]*data[firstVelocity + 2*ie + 1];
      const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
      const CFreal R_gas = K_gas/_m_i[ie];
      const CFreal c_v = c_p - R_gas;
    
      state[endEM + nbSpecies + 2*nbSpecies + ie] = data[firstSpecies + ie]*rho*(c_v*data[firstTemperature + 4*ie] + 0.5*V2);
    } 
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
  

}; // end of class EulerMFMHD2DCons
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD2DConsT_hh
