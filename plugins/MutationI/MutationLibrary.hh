#ifndef COOLFluiD_Physics_Mutation_MutationLibrary_hh
#define COOLFluiD_Physics_Mutation_MutationLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/Fortran.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/LookupTable2D.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutation {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MutationLibrary.
 *
 * @author Tiago Quintino
 * @author Andrea Lani
 * @author Michel Rasquin
 * @author Janos Molnar
 * @author Marco Panesi
 *
 */
class MutationLibrary : public Framework::PhysicalChemicalLibrary {
public:

  enum LambdaAlgo {LAMBDACG=0,LAMBDAD=1};
  enum EtaAlgo    {ETACG=0,ETAD=1};

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Common::LookupTable2D<CFdouble, CFdouble, CFdouble> LkpTable;

  /**
   * Constructor without arguments
   */
  MutationLibrary(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MutationLibrary();
  
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
   * Get the catalycity factor for N
   */
  void getGammaN(CFreal& m_GN);
  
  /**
   * Get the catalycity factor for O
   */
  void getGammaO(CFreal& m_GO);
  
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
    if (_NS == 5) {
      return _Rgas*(1.5 + Y[2]/MOLARMASSP[2] +
		    Y[3]/MOLARMASSP[3] +
		    Y[4]/MOLARMASSP[4]);
    }
    else {
      cf_assert(_NS == 2);
      return _Rgas*(1.5 + Y[1]/MOLARMASSP[1]);
    }
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
      tmpMass +=  Y[is]/MOLARMASSP[is];
    }
    return 1./tmpMass;
  }

  /**
   * Set the constant of gases in J/(Kg*K)
   */
  void setRiGas(RealVector& Ri)
  {
    cf_assert(Ri.size() == static_cast<CFuint>(_NS));
    for(CFint is = 0; is < _NS; ++is) {
      Ri[is] = _Rgas/MOLARMASSP[is];
    }
  }

  /**
   * Set the IDs of the molecules in the mixture
   */
  void setMoleculesIDs(std::vector<CFuint>& v);
  
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
			 CFreal& eta,
			 CFreal& lambdaTrRo, 
			 RealVector& lambdaInt,
			 RealVector& rhoUdiff) {}
 
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
		    RealVector& tempVib,
		    CFdouble& pressure,
		    CFreal& lambdaTrRo,
		    RealVector& lambdaVib);

  /**
   * Calculates the thermal conductivity by conjugate gradient method method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  void lambdaVibNEQN2(CFreal& temperature,
		      RealVector& tempVib,
		      CFdouble& pressure,
		      CFreal& lambdaTrRo,
		      RealVector& lambdaVib);

  /**
   * Calculates the dynamic viscosity, given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble eta(CFdouble& temp, 
	       CFdouble& pressure,
	       CFreal* tempVib)
  {
    if (_etaAlgo == ETACG) return etaCG(temp, pressure, tempVib);
    return etaD(temp, pressure, tempVib);
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
				RealVector* tempVib);

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
      X[i] = x[i];
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
				RealVector& tempVib,
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
   * Sets the species (molar) fractions. This function should be called before getting
   * thermodynamic quantities or transport properties.
   * @param xs the RealVector of the molar fractions of species
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
			     RealVector& tempVib,
			     CFdouble& pressure,
			     CFdouble& rho,
			     const RealVector& ys,
			     bool flagJac,
			     RealVector& omega,
                             RealMatrix& jacobian);

  /**
   * Returns the source term for the vibrational relaxation with VT transfer
   * @param temp the mixture temperature
   * @param tempVib the vibrational temperature
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
   * free electron-electronic conservation equation.
   * The analytical Jacobian matrix of the mass production terms can also be
   * computed for the variables p, u, v,( w,) ys
   * @param pressure the mixture pressure
   * @param temp the mixture temperature
   * @param ys the species mass fractions
   * @param flg_jag the flag to compute the analytical Jacobian matrix
   * @param omega the mass production terms
   * @param jacobian the Jacobian matrix of the mass production terms
   */
  void getSource(CFdouble& temp,
			     RealVector& tempVib,
			     CFdouble& pressure,
			     CFdouble& rho,
			     const RealVector& ys,
			     bool flagJac,
			     RealVector& omega,
                             RealVector& omegav,
                             CFdouble& omegaRad,
                             RealMatrix& jacobian);
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
  virtual void getSourceEE(CFdouble& temp,
			   RealVector& tVec,
			   CFdouble& pressure,
			   CFdouble& rho,
			   const RealVector& ys,
			   bool flagJac,
			   CFdouble& omegaEE) 
  {
    throw Common::NotImplementedException (FromHere(),"MutationLibrary::getSourceEE()");
  }

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
		   RealMatrix& Dij,
		   RealVector& rhoUdiff);
  
  /**
   * Returns the total enthalpies per unit mass of species
   * @param temp the mixture temperature
   * @param pressure the mixture pressure
   */
  void getSpeciesTotEnthalpies(CFdouble& temp,
                               RealVector& tempVib,
			       CFdouble& pressure,
                               RealVector& hsTot,
			       RealVector* hsVib,
			       RealVector* hsEl);
  
protected: // helper function

  typedef CFdouble (MutationLibrary::*ComputeQuantity)
    (CFdouble&, CFdouble&);

  /**
   * Set up all the look up tables
   */
  void setLookUpTables();

  /**
   * Set up the library sequentially
   */
  void setLibrarySequentially();

  /**
   * Compute lookup tables using the ponter to member functions
   * contained in the given vector
   */
  void setTables(std::vector<ComputeQuantity>& varComputeVec);

  /**
   * Calculates the dynamic viscosity by direct method given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble etaD(CFdouble& temp,
		CFdouble& pressure,
		CFreal* tempVib);

  /**
   * Calculates the dynamic viscosity by conjugate gradient method
   * given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble etaCG(CFdouble& temp,
		 CFdouble& pressure,
		 CFreal* tempVib);

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
  
  /**
   * Calculates the density given temperature and pressure.
   * This function is supplied for convinience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param pressure pressure
   */
  CFdouble density(CFdouble& temp,
                   CFdouble& pressure)
  {
    return density(temp,pressure, CFNULL);
  }
  
  
protected: // data to use to interface FORTRAN77

  /// mapping name variable - idx
  Common::CFMap<std::string, size_t> _nameToIdxVar;

  /// table of lookup tables
  LkpTable _lookUpTables;
  
  /// mixture name
  std::string _mixtureName;

  /// reaction name
  std::string _reactionName;

  /// flag telling if to use the look up tables
  bool _useLookUpTable;

  /// flag telling to ignore the electronic energy
  bool _noElectEnergy;

  /// flag telling if to use the look up tables
  std::vector<std::string> _lkpVarNames;

  /// Vector of the impinging fluxes
  RealVector m_mcal;

  /// Min temperature in the table
  CFdouble _Tmin;

  /// Max temperature in the table
  CFdouble _Tmax;

  /// Delta temperature in the table
  CFdouble _deltaT;
  
  /// Catalycity factor for N
  CFdouble _GammaN;
  
  /// Catalycity factor for O
  CFdouble _GammaO;

  /// Min pressure in the table
  CFdouble _pmin;

  /// Max pressure in the table
  CFdouble _pmax;

  /// Delta pressure in the table
  CFdouble _deltaP;

  /// Lagrange-Sonine computation option (polynom order???)
  int _sonine;
  
  /// Where to compute transport properties
  int _imod;

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

  /// flag telling if there are electrons
  int NE;

  /// what is this??? (Number of reactions? : Janos)
  int NREA;

 /// Level at which start the Binaries Diffusion coefficients in WR2 
  int IBINIJ;

  /// lenght of work vector WR1
  int LWR1;

  /// lenght of work vector WR2
  int LWR2;

  /// lenght of work vector WR3
  int LWR3;

  /// lenght of work vector WI
  int LWI;

  /// lenght of work vector WC
  int LWC;

  /// Small disturbance
  CFdouble EPS;

  /// Small disturbance
  CFdouble TOL;
  
  /// Vector of the molar mass
  RealVector m_mmi;

  /// partial derivative of the density in the pressure
  CFdouble DRHODP;

  /// first real valued work vector
  FDOUBLE WR1;

  /// second real valued work vector
  FDOUBLE WR2;

  /// third real valued work vector (for non-equilibrium)
  FDOUBLE WR3;

  /// integer valued work vector
  FINT    WI;

  /// char valued work vector
  FCHAR   WC;

  /// correction functions of the Stefan-Maxwell equations
  FDOUBLE FIJ;

  /// the composition of species
  FDOUBLE X;

  /// the molar composition of species
  FDOUBLE XTOL;

  /// the molar composition initial guess for species
  FDOUBLE Xini;

  /// the nuclear (elemental) mole fractions
  FDOUBLE Xn;

  /// the nuclear (elemental) mass fractions
  FDOUBLE Yn;

  /// the species mass fractions
  FDOUBLE Y;

  /// the mass composition of species initial guess
  FDOUBLE Yini;

  /// elmental heat transfer coefficients
  FDOUBLE LAMBDAEL;

  /// elemental multicomponent difussion coefficients times mixture density
  FDOUBLE ELDIFCOEF;

  /// elemental thermal demixing coefficients times mixture density
  FDOUBLE ELTDIFCOEF;

  /// reaction rates for non-equilibrium
  FDOUBLE OMEGA;

  /// vibrational-translational energy transfer for a given species
  FDOUBLE OMEGAVT;

  /// species molar masses
  FDOUBLE MOLARMASSP;

  /// Driving force: species cell normal concentration gradients
  FDOUBLE DF;

  /// Mass diffusion fluxes of species
  FDOUBLE JDIF;

  /// Mass production terms differentiated by pressure
  FDOUBLE DWDP;

  /// Mass production terms differentiated by temperature
  FDOUBLE DWDT;

  /// Mass production terms differentiated by species mass fractions
  FDOUBLE DWDYI;

  /// total mass enthalpy
  FDOUBLE _HMTOTAL;

  /// translational mass enthalpy
  FDOUBLE _HMTRANS;

  /// electronic mass enthalpy
  FDOUBLE _HMELECT;

  /// rotation mass enthalpy
  FDOUBLE _HMROT;

  /// vibration mass enthalpy
  FDOUBLE _HMVIBR;

  /// formation mass enthalpy
  FDOUBLE _HMFORM;

  /// total enthalpy per mol
  FDOUBLE _HTOTAL;

  /// perturbed total enthalpy per mol
  FDOUBLE _HTOTALP;

  /// translational enthalpy per mol
  FDOUBLE _HTRANS;

  /// perturbed translational enthalpy per mol
  FDOUBLE _HTRANSP;

  /// electronic enthalpy per mol
  FDOUBLE _HELECT;

  /// perturbed electronic enthalpy per mol
  FDOUBLE _HELECTP;

  /// rotation enthalpy per mol
  FDOUBLE _HROT;

  /// perturbed rotation enthalpy per mol
  FDOUBLE _HROTP;

  /// vibration enthalpy per mol
  FDOUBLE _HVIBR;

  /// perturbed vibration enthalpy per mol
  FDOUBLE _HVIBRP;

  /// formation enthalpy per mol
  FDOUBLE _HFORM;

  /// perturbed enthalpy per mol
  FDOUBLE _HFORMP;

  /// total energy per mol
  FDOUBLE _ETOTAL;

  /// translational energy per mol
  FDOUBLE _ETRANS;

  /// electronic energy per mol
  FDOUBLE _EELECT;

  /// rotation energy per mol
  FDOUBLE _EROT;

  /// vibration energy per mol
  FDOUBLE _EVIBR;

  /// formation energy per mol
  FDOUBLE _EFORM;

  /// array of vibrational temperatures
  /// @TODO remove this once the number of Tvibs is equal in coolfluid and
  // in the allocated array inside Mutation
  FDOUBLE _TVARRAY;

  /// array of vibrational energy source terms
  /// @TODO remove this once the number of Tvibs is equal in coolfluid and
  /// in the allocated array inside Mutation
  FDOUBLE _STVIB;

  FDOUBLE _CPE;

  FDOUBLE _CPR;

  FDOUBLE _CPV;

  FDOUBLE _CPINT;

  FDOUBLE _CHIH;

  /// electric charge of species
  FINT _CHARGE;

  /// array of vibrational temperatures PERTURBATIONS
  /// @TODO remove this once the number of Tvibs is equal in coolfluid and
  // in the allocated array inside Mutation
  FDOUBLE _TVIBEPS1;

  FDOUBLE _TVIBEPS;

  FDOUBLE _LAMBDAVIB;

  FDOUBLE _CPVIB;
  
    
  // array for temporary storage
  RealVector _tmp;

  /// Minimum threshold temperature for chemistry
  CFdouble _TminFix;
  
  /// molecules IDs
  std::vector<CFuint> _moleculesIDs;
  
  /// flag indicating molecules IDs
  std::vector<bool> _flagMoleculesIDs;
 
  /// Fator to reduce the stiffness of source terms <1
  CFdouble _factorOmega;


  
}; // end of class MutationLibrary

//////////////////////////////////////////////////////////////////////////////

    } // namespace Mutation

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutation_MutationLibrary_hh
