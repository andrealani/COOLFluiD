#ifndef COOLFluiD_Physics_NitrogenNASA_NitrogenNASALibrary_hh
#define COOLFluiD_Physics_NitrogenNASA_NitrogenNASALibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/Fortran.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NitrogenNASA {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NitrogenNASALibrary.
 *
 * @author Alessandro Munafo'  
 */
class NitrogenNASALibrary : public Framework::PhysicalChemicalLibrary {
public:

  enum LambdaAlgo {LAMBDACG=0,LAMBDAD=1};
  enum EtaAlgo    {ETACG=0,ETAD=1};

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  NitrogenNASALibrary(const std::string& name);

  /**
   * Default destructor
   */
  ~NitrogenNASALibrary();

  /**
   * Setups the path name
   */
  void setLibPathName(const std::string& libPathName)
  {
    _libPath = libPathName;
  }

  /**
   * Setups the mixture name
   */
  void setMixtureName(const std::string& mixtureName)
  {

  }

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Setups the data of the library
   */
  void setup();

  /**
   * Unsetups the data of the library
   */
  void unsetup();
  
  /**
   * Get the molar masses
   */
  void getMolarMasses(RealVector& mm);

  /**
   * Compute and get the electron pressure
   */
  CFdouble electronPressure(CFreal rhoE,
			    CFreal tempE);
  
  /**
   * Get the translational-rotational cv
   * @pre it assumes that the mass fractions
   *      have been already set
   */
  CFdouble getCvTr() const
  {
    CFdouble sumCvTr = 0.0;
    for (CFint i = 0; i < _NS; ++i) {
         sumCvTr += _Y[i]/_MOLARMASSP[i];
    }
    return _Rgas*sumCvTr;
  }

  /**
   * Get the translational-rotational cv
   * @pre it assumes that the mass fractions
   *      have been already set
   */
  CFdouble getMMass() const
  {
    CFdouble tmpMass = 0.0;
    for(CFint is = 0; is < _NS; ++is) {
      tmpMass += _Y[is]/_MOLARMASSP[is];
    }
    return 1./tmpMass;
  }

  /**
   * Set the constant of gases in J/(Kg*K)
   */
  void setRiGas(RealVector& Ri)
  {
    assert(Ri.size() == static_cast<CFuint>(_NS));
    for(CFint is = 0; is < _NS; ++is) {
      Ri[is] = _Rgas/_MOLARMASSP[is];
    }
  }

  /**
   * Set the IDs of the molecules in the mixture
   */
  void setMoleculesIDs(std::vector<CFuint>& v)
  {

  }
  
  /// Set the thermodynamic state (temperature and pressure)
  /// @param species partial densities
  /// @param mixture temperature
  virtual void setState(CFdouble* rhoi, CFdouble* T) {}
  
  /**
   * Calculates the static pressure of the mixture
   * @param rho  density
   * @param temp temperature
   */
  CFdouble pressure(CFdouble& rho,
		    CFdouble& temp,
		    CFreal* tVec);
 
  /**
   * Returns the diffusion flux
   * This function returnm the binary coefficients of Fick
   * But since they are a diagonal matrix it returns as well a scalar
   * It returns as well the diffusive flux 
   * @param temp the mixture temperature
   * @param pressure the mixture pressure

   */ 
  void getDij_fick(RealVector& dx,
		   CFdouble& pressure,
		   CFdouble& temperature,
		   RealMatrix& Dij,
		   RealVector& rhoUdiff); 
  
  /// Computes all transport coefficients for TCNEQ case
  void transportCoeffNEQ(CFreal& temp, 
			 CFdouble& pressure,
			 CFreal* tVec, 
			 RealVector& normConcGradients,
			 RealVector& normTempGradients,
			 CFreal& eta,
			 CFreal& lambdaTrRo, 
			 RealVector& lambdaInt,
			 RealVector& rhoUdiff);
  
  /**
   * Calculates the thermal conductivity by conjugate gradient method method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble lambdaNEQ(CFdouble& temp, CFdouble& pressure);

  /**
   * Calculates the thermal conductivity by conjugate gradient method method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  void lambdaVibNEQ(CFreal& temp,
		    RealVector& tVec,
		    CFdouble& pressure,
		    CFreal& lambdaTrRo,
		    RealVector& lambdaInt);
  
  /**
   * Calculates the dynamic viscosity, given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble eta(CFdouble& temp,
	       CFdouble& pressure,
	       CFreal* tVec)
  {
    if (_etaAlgo == ETACG) return etaCG(temp, pressure, tVec);
    return etaD(temp, pressure, tVec);
  }
  
  /**
   * Calculates the lambda viscosity, given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble lambdaEQ(CFdouble& temp, CFdouble& pressure)
  {
    if (_lambdaAlgo == LAMBDACG) return lambdaCG(temp, pressure);
    return lambdaD(temp, pressure);
  }

  /**
   * Calculates the electrical conductivity given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble sigma(CFdouble& temp,
     CFdouble& pressure,
     CFreal* tVec);
  // Vatsalya
 virtual CFdouble sigma_debug(CFdouble& temp, CFdouble& pressure, CFreal* tVec, CFuint elem_no){ 
    return 0;
  }
  /**
   * Calculates the specific heat ratio and the speed of sound in
   * thermal equilibrium.
   * @param temp temperature
   * @param pressure pressure
   * @param rho density
   * @param gamma specific heat ratio
   * @param soundSpeed speed
   */
  void gammaAndSoundSpeed(CFdouble& temp,
        CFdouble& pressure,
        CFdouble& rho,
        CFdouble& gamma,
        CFdouble& soundSpeed);

  /**
   * Calculates the specific heat ratio and the speed of sound in
   * thermal equilibrium.
   * @param temp temperature
   * @param pressure pressure
   * @param rho density
   * @param gamma specific heat ratio
   * @param soundSpeed speed
   */
  void frozenGammaAndSoundSpeed(CFdouble& temp,
        CFdouble& pressure,
        CFdouble& rho,
        CFdouble& gamma,
        CFdouble& soundSpeed,
        RealVector* tVec);

  /**
   * Calculates the composition given temperature and pressure.
   * @param temp temperature
   * @param pressure pressure
   * @param x composition array (one component for each species)
   * @pre this function ALWAYS computes the composition
   *      (no lookup table based results...)
   */
  void setComposition(CFdouble& temp,
          CFdouble& pressure,
          RealVector* x);

  /**
   * Reset the composition
   */
  void resetComposition(const RealVector& x)
  {
    for (CFint i = 0; i < _NS; ++i) {
      _X[i] = x[i];
    }
  }

  /**
   * Calculates the density, the enthalpy and the internal energy
   * This function is supplied for convenience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param tempV vibrational temperature
   * @param pressure pressure
   * @param dhe array with density, enthalpy, energy. vibrational energies
   * @param storeExtraData  flag telling if extra data have to be stored
   */
  void setDensityEnthalpyEnergy(CFdouble& temp,
        RealVector& tVec,
        CFdouble& pressure,
        RealVector& dhe,
        bool storeExtraData = false);

  /**
   * Calculates the density, the enthalpy and the internal energy
   * This function is supplied for convenience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param pressure pressure
   */
  void setDensityEnthalpyEnergy(CFdouble& temp,
        CFdouble& pressure,
        RealVector& dhe);


  /**
   * Calculates the density given temperature and pressure.
   * This function is supplied for convinience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble density(CFdouble& temp,
       CFdouble& pressure,
       CFreal* tVec);

  /**
   * Calculates the internal energy at given temperature
   * and pressure.
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble energy(CFdouble& temp,
      CFdouble& pressure);

  /**
   * Calculates the energy in LTE conditions
   * at given temperature and pressure.
   * @param temp      temperature
   * @param pressure  pressure
   * @param intEnergy computed internal energy
   */
  CFdouble enthalpy(CFdouble& temp,
        CFdouble& pressure);

  /**
   * Calculates the speed of sound in
   * thermal equilibrium.
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble soundSpeed(CFdouble& temp,
          CFdouble& pressure);

  /**
   * Sets the mole fractions of elements (nucleons) Xn for the NitrogenNASA
   * environment (for variable elemental composition), this function
   * should be called before getting properties related to the elemental
   * fractions in case of LTE with Demixing
   * @param yn the RealVector of the mass fractions of elements
   */
  void setElemFractions(const RealVector& yn);

  /**
   * Sets the mole fractions of elements (nucleons) Xn for the NitrogenNASA
   * environment (for variable elemental composition) starting from the given
   * species mass fractions
   * @param yn the RealVector of the mass fractions of species
   */
  void setElementXFromSpeciesY(const RealVector& ys);

  /**
   * Sets the species molar fractions. This function should be called before getting
   * thermodynamic quantities or transport properties.
   * @param ys the RealVector of the mass fractions of species
   */
  void setSpeciesFractions(const RealVector& ys);

  /**
   * Sets the species molar fractions. This function should be called before getting
   * thermodynamic quantities or transport properties.
   * @param ys the RealVector of the mass fractions of species
   */
  void setSpeciesMolarFractions(const RealVector& xs);

  /**
   * Sets the electron fractions in the mass composition according to charge
   * neutrality.
   * @param ys the RealVector of the mass fractions of species
   */
  void setElectronFraction(RealVector& ys);

  /**
   * Gets the species (molar) fractions.
   * @param ys the RealVector of the mass fractions of species
   * @param xs the RealVector of the molar fractions of species
   */
  void getSpeciesMolarFractions(const RealVector& ys, RealVector& xs);

  /**
   * Gets the species (mass) fractions.
   * @param xs the RealVector of the molar fractions of species
   * @param ys the RealVector of the mass fractions of species
   */
  void getSpeciesMassFractions(const RealVector& xs, RealVector& ys);

  /**
   * Gets the species mass fractions.
   * @param ys the RealVector of the mass fractions of species
   */
  void getSpeciesMassFractions(RealVector& ys);

  /**
   * Returns the transport coefficients for the diffusive fluxes used in
   * computations with LTE and Demixing
   * @param temp the mixture temperature
   * @param pressure the mixture pressure
   * @param lambda the Butler and Brokaw thermal reactive conductivity
   * @param lambdacor the demixing thermal reactive conductivity
   * @param lambdael the elmental heat transfer coefficients
   * @param eldifcoef elemental multicomponent difussion coefficients
   *        times mixture density
   * @param eltdifcoef the elemental thermal demixing coefficients
   *        times mixture density
   */
  void getTransportCoefs(CFdouble& temp,
			 CFdouble& pressure,
			 CFdouble& lambda,
			 CFdouble& lambdacor,
			 RealVector& lambdael,
			 RealMatrix& eldifcoef,
			 RealVector& eltdifcoef);
  
  /**
   * Returns the mass production/destruction terms [kg m^-3 s^-1] in CHEMICAL
   * NONEQUILIBRIUM based on Arrhenius's formula.
   * The analytical Jacobian matrix of the mass production terms can also be
   * computed for the variables p, u, v,( w,) ys
   * @param pressure the mixture pressure
   * @param temp the mixture temperature
   * @param ys the species mass fractions
   * @param flg_jag the flag to compute the analytical Jacobian matrix
   * @param omega the mass production terms
   * @param jacobian the Jacobian matrix of the mass production terms
   */
  void getMassProductionTerm(CFdouble& temp,
			     RealVector& tVec,
			     CFdouble& pressure,
			     CFdouble& rho,
			     const RealVector& ys,
			     bool flagJac,
			     RealVector& omega,
                             RealMatrix& jacobian);
  
  /**
   * Returns the source term for the vibrational relaxation with VT transfer
   * @param temp the mixture temperature
   * @param tVec the vibrational temperature
   * @param pressure the mixture pressure
   * @param rho the mixture density
   * @param omegav the source term
   */
  void getSourceTermVT(CFdouble& temp,
		       RealVector& tVec,
		       CFdouble& pressure,
		       CFdouble& rho,
		       RealVector& omegav,
		       CFdouble& omegaRad);
  
  /**
   * Returns the source terms species continuity equations, 
   * vibrational energy conservation equation and 
   * free electron-electronic conservation equation
   * @param temp the mixture temperature
   * @param tVec the vibrational temperature
   * @param pressure the mixture pressure
   * @param ys the species mass fractions
   * @param flg_jag the flag to compute the analytical Jacobian matrix
   * @param rho the mixture density
   * @param omegaEE the source term
   */
  void getSourceEE(CFdouble& temp,
		   RealVector& tVec,
		   CFdouble& pressure,
		   CFdouble& rho,
		   const RealVector& ys,
		   bool flagJac,
		   CFdouble& omegaEE);
  
   /**
   * Returns the source terms species continuity equations, 
   * vibrational energy conservation equation and 
   * free electron-electronic conservation equation.
   * The analytical Jacobian matrix of the mass production terms can also be
   * computed for the variables p, u, v,( w,) ys 
   * @param temp the mixture temperature
   * @param tVec the vibrational temperature
   * @param pressure the mixture pressure
   * @param ys the species mass fractions
   * @param flg_jag the flag to compute the analytical Jacobian matrix
   * @param rho the mixture density
   * @param omega the mass production terms
   * @param omegav the source term
   * @param jacobian the Jacobian matrix of the mass production terms
   */
   void getSource(CFdouble& temp,
		 RealVector& tVec,
		 CFdouble& pressure,
		 CFdouble& rho,
		 const RealVector& ys,
		 bool flagJac,
		 RealVector& omega,
                 RealVector& omegav,
                 CFdouble& omegaRad,
                 RealMatrix& jacobian);

  /**
   * Returns the diffusion velocities of species multiplied by the species
   * densities for nonequilibrium computations
   * @param temp the mixture temperature
   * @param pressure the mixture pressure
   * @param normConcGradients the cell normal gradients of species mass fractions
   */
  void getRhoUdiff(CFdouble& temp,
                   CFdouble& pressure,
                   RealVector& normConcGradients,
		   RealVector& normTempGradients,
		   CFreal* tVec,
                   RealVector& rhoUdiff,
		   bool fast);

  /**
   * Returns the total enthalpies per unit mass of species
   * @param temp the mixture temperature
   * @param pressure the mixture pressure
   */
  void getSpeciesTotEnthalpies(CFdouble& temp,
                               RealVector& tVec,
			       CFdouble& pressure,
                               RealVector& hsTot,
			       RealVector* hsVib,
			       RealVector* hsEl);
  
private: // helper function

   /**
    * Copy data files
    */
    void copyDataFiles();

   /**
    * Delete data files
    */
    void deleteDataFiles();

  /**
   * Test and write thermodynamic and transport properties
   */
  void testAndWriteProperties();

  /**
   * Test and write the source terms
   */
  void testAndWriteSourceTerms();

  /**
   * Temperature of free electrons
   */
  CFdouble getTe(CFdouble temp, CFreal* tVec)
  {
    return 0;
  }
  
  /**
   * Set the working temperature vector
   */
  void setTVarray(RealVector& tVec)
  {

  }

  /**
   * Set the working temperature vector
   */
  void setTVec(CFdouble& tr, RealVector& tVec)
  {

  }
    
  /**
   * Set the working temperature vector
   */
  void setDefaultTVarray(CFreal& temp)
  {
    
  }

  /**
   * Set up the library sequentially
   */
  void setLibrarySequentially();

  /**
   * Calculates the dynamic viscosity by direct method given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble etaD(CFdouble& temp,
    CFdouble& pressure,
    CFreal* tVec);

  /**
   * Calculates the dynamic viscosity by conjugate gradient method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble etaCG(CFdouble& temp,
     CFdouble& pressure,
     CFreal* tVec);

  /**
   * Calculates the thermal conductivity by direct method method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble lambdaD(CFdouble& temp,
       CFdouble& pressure);

  /**
   * Calculates the thermal conductivity by conjugate gradient method method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble lambdaCG(CFdouble& temp,
        CFdouble& pressure);

private: // data to use to interface FORTRAN77

  /// library path
  std::string _libPath;

  /// flag for the model to be used (vibrational specific or multi-temperature) 
  std::string _libPhysModel;

  /// flag for the vibrational specific model to be used 
  std::string _libFlagVibSTS;

  /// Where to compute transport properties
  CFint _imod;

  /// Small value for mole fractions to limit the composition with COMPOTOL2
  CFdouble _Xlim;

  /// algorithm to compute thermal conductivity
  std::string _lambdaAlgoStr;

  /// algorithm to compute dynamic viscosity
  std::string _etaAlgoStr;

  /// algorithm to compute thermal conductivity
  LambdaAlgo _lambdaAlgo;

  /// algorithm to compute dynamic viscosity
  EtaAlgo    _etaAlgo;

  /// Small disturbance
  CFdouble _EPS;

  /// Small disturbance
  CFdouble _TOL;

  /// Species densities
  FDOUBLE _Rhoi;

  /// the composition of species
  FDOUBLE _X;

  /// the molar composition of species
  FDOUBLE _XTOL;

  /// the molar composition initial guess for species
  FDOUBLE _Xini;

  /// the nuclear (elemental) mole fractions
  FDOUBLE _Xn;

  /// the nuclear (elemental) mass fractions
  FDOUBLE _Yn;

  /// the species mass fractions
  FDOUBLE _Y;

  /// the mass composition of species initial guess
  FDOUBLE _Yini;

  /// reaction rates for non-equilibrium
  FDOUBLE _OMEGA;

   /// species molar masses
  FDOUBLE _MOLARMASSP;

  /// Mass production terms differentiated by species mass fractions
  FDOUBLE _OMEGAJACOB;
  
  /// jacobian of driving forces for vibrational and electronic energies
  FDOUBLE _OMEGAVTJACOB;
    
  /// total enthalpy per mass
  FDOUBLE _HTOTALUM;

  /// translational enthalpy per mass
  FDOUBLE _HTRANSUM;

  /// electronic enthalpy per mass
  FDOUBLE _HELECTUM;

  /// rotation enthalpy per mass
  FDOUBLE _HROTUM;

  /// vibration enthalpy per mass
  FDOUBLE _HVIBRUM;

  /// formation enthalpy per mass
  FDOUBLE _HFORMUM;

  /// array of vibrational energy source terms
  FDOUBLE _OMEGAVIB;

  /// Constant volume frozen specific vector
  FDOUBLE _CV;

  /// Constant volume frozen specific vector
  FDOUBLE _CP;

  /// Minimum threshold temperature for chemistry
  CFdouble _TminFix;
 
}; // end of class NitrogenNASA

//////////////////////////////////////////////////////////////////////////////

    } // namespace NitrogenNASA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NitrogenNASA_NitrogenNASA_hh
