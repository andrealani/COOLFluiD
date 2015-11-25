#ifndef COOLFluiD_Physics_Mutation2OLD_MutationLibrary2OLD_hh
#define COOLFluiD_Physics_Mutation2OLD_MutationLibrary2OLD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/Fortran.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutation2OLD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MutationLibrary2OLD.
 *
 * @author Tiago Quintino
 * @author Andrea Lani
 * @author Michel Rasquin
 * @author Janos Molnar
 * @author Marco Panesi
 *
 */
class MutationLibrary2OLD : public Framework::PhysicalChemicalLibrary {
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
  MutationLibrary2OLD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MutationLibrary2OLD();
  
  /**
   * Get the catalycity factor for N
   **/
  void getGammaN(CFreal& m_GN);

  /**
   * Get the catalycity factor for O
   **/
  void getGammaO(CFreal& m_GO);
  
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
  virtual void setup();
  
  /**
   * Unsetups the data of the library
   */
  virtual void unsetup();

  /**
   * Get the molar masses
   */
  virtual void getMolarMasses(RealVector& mm);

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
    CFdouble sumCvTr = 1.5;
    for (CFint i = 0; i < _NS; ++i) {
      if (_WI[_IATOMI+i-1] > 1) {
	sumCvTr += _Y[i]/_MOLARMASSP[i];
      }
    }
    return _Rgas*sumCvTr;
  }

  /**
   * Get the mixture molar mass 
   * @pre it assumes that the mass fractionshave been already set
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
    v.reserve(_NV);
    for (CFint i = 0; i < _NS; ++i) {
      if (_WI[_IATOMI+i-1] > 1) {
	v.push_back(i);
      }
    }
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
  virtual CFdouble eta(CFdouble& temp,
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
  virtual CFdouble lambdaEQ(CFdouble& temp, CFdouble& pressure)
  {
    if (_lambdaAlgo == LAMBDACG) return lambdaCG(temp, pressure);
    return lambdaD(temp, pressure);
  }

  /**
   * Calculates the electrical conductivity given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  virtual CFdouble sigma(CFdouble& temp,
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
  virtual void gammaAndSoundSpeed(CFdouble& temp,
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
  virtual void frozenGammaAndSoundSpeed(CFdouble& temp,
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
 virtual void setComposition(CFdouble& temp,
			     CFdouble& pressure,
			     RealVector* x);
  
  /**
   * Reset the composition
   */
  virtual void resetComposition(const RealVector& x)
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
  virtual void setDensityEnthalpyEnergy(CFdouble& temp,
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
  virtual void setDensityEnthalpyEnergy(CFdouble& temp,
					CFdouble& pressure,
					RealVector& dhe);
  

  /**
   * Calculates the density given temperature and pressure.
   * This function is supplied for convinience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param pressure pressure
   */
  virtual CFdouble density(CFdouble& temp,
			   CFdouble& pressure,
			   CFreal* tVec);
  
  /**
   * Calculates the internal energy at given temperature
   * and pressure.
   * @param temp temperature
   * @param pressure pressure
   */
  virtual CFdouble energy(CFdouble& temp,
			  CFdouble& pressure);
  
  /**
   * Calculates the energy in LTE conditions
   * at given temperature and pressure.
   * @param temp      temperature
   * @param pressure  pressure
   * @param intEnergy computed internal energy
   */
  virtual  CFdouble enthalpy(CFdouble& temp,
			     CFdouble& pressure);
  
  /**
   * Calculates the speed of sound in
   * thermal equilibrium.
   * @param temp temperature
   * @param pressure pressure
   */
  virtual CFdouble soundSpeed(CFdouble& temp,
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
   * Returns the source terms species continuity equations, 
   * vibrational energy conservation equation and 
   * free electron-electronic conservation equation
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
  
  /// Computes all transport coefficients for TCNEQ case
  void transportCoeffNEQ(CFreal& temp, 
			 CFdouble& pressure,
			 CFreal* tVec, 
			 RealVector& normConcGradients,
			 CFreal& eta,
			 CFreal& lambdaTrRo, 
			 RealVector& lambdaInt,
			 RealVector& rhoUdiff);
  
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

protected: // helper function
  
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
   * Set the working temperature vector
   */
  void setTVarray(RealVector& tVec)
  {
    using namespace std;
    
    CFuint ic = 0;
    for (CFint i = 0; i < _NS; ++i) {
      const CFuint nvibMode = _WI[_IVIBI+i-1];
      for (CFuint j = 0; j < nvibMode; ++j, ++ic) {
	_TVARRAY[ic] = tVec[_WI[_IVIBTEMPI+i-1]-2];
      }
    }
  }

  /**
   * Set the working temperature vector
   */
  void setTVec(CFdouble& tr, RealVector& tVec)
  {
    _TVEC[0]= tr;
    for (CFint i = 0; i < _nbTvib; ++i) {
      _TVEC[1+i] = tVec[i];
    }
    _TVEC[1+_nbTvib] = (_nbTe == 0) ? tVec[0] : tVec[tVec.size()-1];
  }

  /**
   * Set the working temperature vector
   */
  void setDefaultTVarray(CFreal& temp)
  {
    for (CFint i = 0; i < _NV; ++i) {
      _TVARRAY[i] = temp;
    }
  }

  /**
   * Candler's CV model
   */
  void cvCandler();

  /**
   * Gnoffo's CV model
   */
  void cvGnoffo();

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

protected: // data to use to interface FORTRAN77
  
  /// mixture name
  std::string _mixtureName;

  /// reaction name
  std::string _reactionName;

  /// transfer file name
  std::string _transfName;
  
  /// thermo dir
  std::string _thermoDir;
  
  /// flag telling to ignore the electronic energy
  bool _noElectEnergy;

  /// include the electronic energy
  bool _includeElectronicEnergy;

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
  int _NE;

  /// flag telling the number of molecules
  int _NV;

  /// Number of reactions
  int _NREA;
  
  /// Level at which start the Binaries Diffusion coefficients in WR2
  int _IBINIJ;
  
  /// lenght of work vector WR1
  int _LWR1;

  /// lenght of work vector WR2
  int _LWR2;

  /// lenght of work vector WR3
  int _LWR3;

  /// lenght of work vector WR4
  int _LWR4;

  /// lenght of work vector WI
  int _LWI;

  /// lenght of work vector WC
  int _LWC;

  ///  Index of vibrational temperatures per species
  int _IVIBTEMPI;

  ///  Index of species per each equation (vibrational Temperature)
  int _IVIBSPEI;

  ///  NUMBER of species per each equation (vibrational Temperature) I.E.(TVIB(1) == N2+NOp etc)
  int _INVIBSPEI;

  ///  NUMBER of species per each equation (vibrational Temperature) I.E.(TVIB(1) == N2+NOp etc)
  int _IVIBI;

  ///  ATOMS per species
  int _IATOMI;

  ///  index for dissociation
  int _IDIS;

  /// Small disturbance
  CFdouble _EPS;

  /// Small disturbance
  CFdouble _TOL;

  /// partial derivative of the density in the pressure
  CFdouble _DRHODP;

  /// first real valued work vector
  FDOUBLE _WR1;

  /// second real valued work vector
  FDOUBLE _WR2;

  /// third real valued work vector (for non-equilibrium)
  FDOUBLE _WR3;

  /// fourth real valued work vector (for non-equilibrium)
  FDOUBLE _WR4;

  /// integer valued work vector
  FINT _WI;

  /// char valued work vector
  FCHAR _WC;

  /// correction functions of the Stefan-Maxwell equations
  FDOUBLE _FIJ;

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

  /// elmental heat transfer coefficients
  FDOUBLE _LAMBDAEL;

  /// elemental multicomponent difussion coefficients times mixture density
  FDOUBLE _ELDIFCOEF;

  /// elemental thermal demixing coefficients times mixture density
  FDOUBLE _ELTDIFCOEF;

  /// reaction rates for non-equilibrium
  FDOUBLE _OMEGA;

  /// reaction rates for ionization non-equilibrium
  CFdouble _OMEGAI;

  /// species molar masses
  FDOUBLE _MOLARMASSP;

  /// Driving force: species cell normal concentration gradients
  FDOUBLE _DF;

  /// Mass diffusion fluxes of species
  FDOUBLE _JDIF;

  /// Mass production terms differentiated by species mass fractions
  FDOUBLE _OMEGAJACOB;
  
  /// jacobian of driving forces for vibrational and electronic energies
  FDOUBLE _OMEGAVTJACOB;
  
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

  /// total mass enthalpy
  FDOUBLE _HT;

  /// translational mass enthalpy
  FDOUBLE _HTT;

  /// electronic mass enthalpy
  FDOUBLE _HTE;

  /// rotation mass enthalpy
  FDOUBLE _HTR;

  /// vibration mass enthalpy
  FDOUBLE _HTV;

  /// formation mass enthalpy
  FDOUBLE _HTF;

  /// total mass enthalpy
  FDOUBLE _HE;

  /// translational mass enthalpy
  FDOUBLE _HET;

  /// electronic mass enthalpy
  FDOUBLE _HEE;

  /// rotation mass enthalpy
  FDOUBLE _HER;

  /// vibration mass enthalpy
  FDOUBLE _HEV;

  /// formation mass enthalpy
  FDOUBLE _HEF;

  /// total enthalpy per mol
  FDOUBLE _HTOTAL;

  /// translational enthalpy per mol
  FDOUBLE _HTRANS;

  /// electronic enthalpy per mol
  FDOUBLE _HELECT;

  /// rotation enthalpy per mol
  FDOUBLE _HROT;

  /// vibration enthalpy per mol
  FDOUBLE _HVIBR;

  /// formation enthalpy per mol
  FDOUBLE _HFORM;
  
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

  /// perturbed total enthalpy per mol
  FDOUBLE _HTOTALP;

  /// perturbed translational enthalpy per mol
  FDOUBLE _HTRANSP;

  /// perturbed electronic enthalpy per mol
  FDOUBLE _HELECTP;

  /// perturbed rotation enthalpy per mol
  FDOUBLE _HROTP;

  /// perturbed vibration enthalpy per mol
  FDOUBLE _HVIBRP;

  /// perturbed enthalpy per mol
  FDOUBLE _HFORMP;

  /// internal enthalpy per mol
  FDOUBLE _HINT;

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

  FDOUBLE _CPE;

  FDOUBLE _CPR;

  FDOUBLE _CPV;

  FDOUBLE _CPINT;

  FDOUBLE _CPVIB;

  FDOUBLE _CHIH;

  /// electric charge of species
  FINT _CHARGE;

  /// array of vibrational temperatures
  FDOUBLE _TVARRAY;

  /// array of vibrational energy source terms
  FDOUBLE _STVIB;

  /// array of vibrational temperatures PERTURBATIONS
  FDOUBLE _TVIBEPS1;

  FDOUBLE _TVIBEPS;

  FDOUBLE _LAMBDAVIB;

  /// SOURCE TERM for chemistry for Thermal non-equilibrium
  FDOUBLE _OMEGACV;

  /// vibrational-translational energy transfer for a given species
  FDOUBLE _OMEGAVE;

  /// vibrational-vibrational transfer for a given species
  FDOUBLE _OMEGAVV;

  FDOUBLE _TVEC;

  FDOUBLE _TT;

  FDOUBLE _TE;

  /// Minimum threshold temperature for chemistry
  CFdouble _TminFix;

  /// Catalycity factor for N
  CFdouble _GammaN;
  
  /// Catalycity factor for O
  CFdouble _GammaO;

  /// escape factor
  CFdouble _escape;
  
  /// iimol factor
  int _iimol;
  
  /// model for the CV coupling
  CFuint _cvModel;

  /// preferential dissociation factors (size == number of molecules)
  std::vector<CFreal> _prefDissFactor;

  /// molecule IDs
  std::vector<CFuint> _molecIDs;
 
  /// flag indicating molecules IDs
  std::vector<bool> _flagMoleculesIDs;
  
  /// IDs of the species whose electronic energy has to be considered
  std::vector<CFuint> _boltzmannIDs;

  /// ID to identify the molecule that gets all Tv contributions
  CFuint _molecTvID;

  /// Factor to reduce the stiffness of source terms <1
  CFdouble _factorOmega;
  
  /// add Ramshaw correction to fick law
  bool _useRamshaw;
  
}; // end of class MutationLibrary2OLD

//////////////////////////////////////////////////////////////////////////////

    } // namespace Mutation2OLD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutation2OLD_MutationLibrary2OLD_hh
