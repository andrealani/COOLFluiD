#ifndef COOLFluiD_Physics_Mutationpp_MutationLibrarypp_hh
#define COOLFluiD_Physics_Mutationpp_MutationLibrarypp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/Fortran.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include <mutation++.h>
#include "Common/LookupTable2D.hh" //@modif_LkT

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutationpp {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MutationLibrarypp.
 *
 * @author Andrea Lani
 * @author James B. Scoggins
 *
 */
class MutationLibrarypp : public Framework::PhysicalChemicalLibrary {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  /**
   * Creation of an alias for the class that defines a table
   * @modif_LkT
   */
  typedef Common::LookupTable2D<CFdouble, CFdouble, CFdouble> LkpTable; 
  /**
   * Constructor without arguments
   */
  MutationLibrarypp(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MutationLibrarypp();
  
  /**
   * Setups the mixture name
   */
  void setMixtureName(const std::string& mixtureName)
  {
    _mixtureName = mixtureName;
  }

  /**
   * Configures this configurable object.
   */
  void configure ( Config::ConfigArgs& args );

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
  
  /// Set the thermodynamic state (temperature and pressure)
  /// @param species partial densities
  /// @param mixture temperature
  void setState(CFdouble* rhoi, CFdouble* T) 
  {
    for (CFuint i = 0; i < _NS; ++i) {
      m_rhoiv[i] = std::max(_minRhoi, rhoi[i]);
    }
    CFLog(DEBUG_MAX, "MutationLibrarypp::setState() => rhoiv = " << m_rhoiv << ", T = " << *T << "\n"); 
    
    // this needs to be fixed for 2-temperatures
    for (CFuint i = 0; i < m_Tstate.size(); ++i) {
      m_Tstate[i] = std::max(T[i], _minT);
    }
    m_gasMixture->setState(&m_rhoiv[0], &m_Tstate[0], 1);
  }   
  
  /**
   * Compute and get the electron pressure
   */
  CFdouble electronPressure(CFreal rhoE, CFreal tempE);
  
   /**
    * Get the translational-rotational cv
    * @pre it assumes that the mass fractions
    *      have been already set
    */
    CFdouble getCvTr() const
   {
     throw Common::NotImplementedException(FromHere(),"MutationLibrarypp::getCvTr()");
     return 0.;
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
      tmpMass += m_y[is]/m_molarmassp[is];
    }
    return 1./tmpMass;
  }
  
   /**
    * Set the constant of gases in J/(Kg*K)
    */
  void setRiGas(RealVector& Ri)
  {
    assert(Ri.size() == (CFuint)_NS);
    cf_assert(m_molarmassp.size() == (CFuint)_NS);
    for(CFint is = 0; is < _NS; ++is) {
      Ri[is] = Mutation::RU/m_molarmassp[is];
    }
  }

   /**
    * Set the IDs of the molecules in the mixture
    */
   void setMoleculesIDs(std::vector<CFuint>& v)
  {
    // AL: this is buggy, it always put electron first
    const CFuint nbSpecies = getNbSpecies();
    for (CFuint i = 0; i < nbSpecies; ++i) {
      if (m_gasMixture->species(i).type() == Mutation::Thermodynamics::MOLECULE) {
	CFLog(DEBUG_MAX, "MutationLibrarypp::setMoleculesIDs() => species ["<< i << "] is " 
	      << m_gasMixture->species(i).name() << "\n");
	v.push_back(i);
      }
    }
  }
  
  /**
   * Calculates the static pressure of the mixture
   * @param rho  density
   * @param temp temperature
   */
  CFdouble pressure(CFdouble& rho,
		    CFdouble& temp,
		    CFreal* tVec);
  
  
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
  CFdouble eta(CFdouble& temp, CFdouble& pressure, CFreal* tVec)
  {
    //AL: here make sure that if Te is present, viscosity() uses that one 
    CFreal mu = m_gasMixture->viscosity();
    CFLog(DEBUG_MAX, "Mutation::eta() => mu = " << mu << "\n");
    return mu;
  }
  
  /**
   * Calculates the lambda viscosity, given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
    CFdouble lambdaEQ(CFdouble& temp, CFdouble& pressure)
    {
      return m_gasMixture->equilibriumThermalConductivity();
    }
  
  /**
   * Calculates the electrical conductivity given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble sigma(CFdouble& temp,
		 CFdouble& pressure,
		 CFreal* tVec);
  
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
       m_x[i] = x[i];
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
//@modif_LkT
  CFdouble density(CFdouble& temp,
                   CFdouble& pressure)
  {
    return density(temp,pressure, CFNULL);
  }
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
   * Sets the mole fractions of elements (nucleons) Xn for the Mutation
   * environment (for variable elemental composition), this function
   * should be called before getting properties related to the elemental
   * fractions in case of LTE with Demixing
   * @param yn the RealVector of the mass fractions of elements
   */
  void setElemFractions(const RealVector& yn);
  
  /**
   * Sets the mole fractions of elements (nucleons) Xn for the Mutation
   * environment (for variable elemental composition) starting from the given
   * species mass fractions
   * @param yn the RealVector of the mass fractions of species
   */
  void setElementXFromSpeciesY(const RealVector& ys);
  
  /**
   * Sets the species (molar) fractions. This function should be called before getting
   * thermodynamic quantities or transport properties.
   * @param ys the RealVector of the mass fractions of species
   */
  void setSpeciesFractions(const RealVector& ys);
  
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
		    CFreal& Diff_coeff,
		    RealMatrix& Dij,
		    RealVector& rhoUdiff);
  /**
   ** Get the catalycity factor for N
   */
  void getGammaN(CFreal& m_GN);
  
  /**
   * Get the catalycity factor for O
   */
   void getGammaO(CFreal& m_GO);
  
   /**
    * Sets the species (molar) fractions. This function should be called before getting
    * thermodynamic quantities or transport properties.
    * @param xs the RealVector of the molar fractions of species
    */
   void setSpeciesMolarFractions(const RealVector& xs);

  
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
    
  /// enumerator for the state model type 
  enum StateModelType {LTE=0, CNEQ=1, TCNEQ=2};

  /* ========================
     START section @modif_LkT 
  */
protected: // helper function
  typedef CFdouble (MutationLibrarypp::*ComputeQuantity)
    (CFdouble&, CFdouble&);


  /**
   * Set up all the look up tables 
   */
  void setLookUpTables();

	/**
	 * Set up the library sequentially 
	 */
	//void setLibrarySequentially();

  /**
   * Compute lookup tables using the ponter to member functions
   * contained in the given vector
   */
  void setTables(std::vector<ComputeQuantity>& varComputeVec);
  
protected: //variables

  /// flag telling if to use the look up tables
  bool _useLookUpTable;

  /// flag telling to ignore the electronic energy
  bool _noElectEnergy;

  /// mapping name variable - idx
  Common::CFMap<std::string, size_t> _nameToIdxVar;

  /// table of lookup tables
  LkpTable _lookUpTables;

  /// gas mixture pointer
  std::auto_ptr<Mutation::Mixture> m_gasMixture;
  
  /// gas mixture pointer for equilibrium
  Mutation::Mixture* m_gasMixtureEquil; 
  
  /// state model type enumerator
  MutationLibrarypp::StateModelType m_smType;

  /// local number of Tv (to be removed)
  CFuint _nbTvibLocal;
  
  /// flag telling if to use the look up tables
  std::vector<std::string> _lkpVarNames;

  /// Vector of the impinging fluxes
  RealVector m_mcal;

  /// Vector of the molar mass
  RealVector m_mmi;

  /// mass fractions
  RealVector m_y;
  
  /// molar fractions
  RealVector m_x;
  
  /// the nuclear (elemental) mass fractions
  RealVector m_yn;
  
  /// the nuclear (elemental) molar fractions
  RealVector m_xn;
  
  /// species molar masses
  RealVector m_molarmassp;
  
  /// stores the charge of each species
  RealVector m_charge;
  
  /// stores the modified driving forces for the Stefan-Maxwell system solution
  RealVector m_df;
  
  /// backup of the partial densities
  RealVector m_rhoivBkp;
  
  /// partial densities
  RealVector m_rhoiv;
  
  /// enthalpies
  RealVector m_ht;

  /// enthalpies
  RealVector m_hr;  

  /// enthalpies
  RealVector m_hf; 

  /// state temperatures
  RealVector m_Tstate;

  /// mixture name
  std::string _mixtureName;
    
  /// state model name
  std::string _stateModelName;

  /// minimum partial density
  CFreal _minRhoi;
  
  /// minimum temperature
  CFreal _minT;

  /// Min temperature in the table
  CFdouble _Tmin;

  /// Max temperature in the table
  CFdouble _Tmax;

  /// Minimum threshold temperature for chemistry
  CFdouble _TminFix;

  /// Delta temperature in the table
  CFdouble _deltaT;

  ///Logarithmic scale for the pressure
  bool _pLogScale;

  /// Min pressure in the table
  CFdouble _pmin;

  /// Max pressure in the table
  CFdouble _pmax;

  /// Delta pressure in the table
  CFdouble _deltaP;

  /// Small disturbance
  CFdouble EPS;

}; // end of class MutationLibrarypp
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace Mutationpp

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutationpp_MutationLibrarypp_hh
