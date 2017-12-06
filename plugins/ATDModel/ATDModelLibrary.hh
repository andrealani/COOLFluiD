#ifndef COOLFluiD_Physics_ATDModel_ATDModelLibrary_hh
#define COOLFluiD_Physics_ATDModel_ATDModelLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalChemicalLibrary.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ATDModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a ATDModelLibrary.
 * For the order of the species in the mixture look *mix at ../Mutation2.0I/data/mixture/
 *
 * @author Francisco Castro Rego
 * @author Jesus Garicano Mena
 * @author Andrea Lani
 *
 */
class ATDModelLibrary : public Framework::PhysicalChemicalLibrary {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  ATDModelLibrary(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ATDModelLibrary();
  
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
   * @pre mm.size() = getNbSpecies()
   */
  virtual void getMolarMasses(RealVector& mm);

  /**
   * Compute and get the electron pressure
   * @param rhoE  electron density
   * @param tempE electron temperature
   */
  virtual CFdouble electronPressure(CFreal rhoE, CFreal tempE);

  /**
   * Get the translational-rotational cv
   * @pre it assumes that the mass fractions have been already set
   */
  virtual CFdouble getCvTr() const;

  /**
   * Get the translational-rotational cv
   * @pre it assumes that the mass fractions
   *      have been already set
   */
  virtual CFdouble getMMass() const;

  /**
   * Set the constant of gases in J/(Kg*K) for each species
   */
  virtual void setRiGas(RealVector& Ri);

  /**
   * Set the IDs of the molecules in the mixture
   */
  virtual void setMoleculesIDs(std::vector<CFuint>& v);
  
  /// Set the thermodynamic state (temperature and pressure)
  /// @param species partial densities
  /// @param mixture temperature
  virtual void setState(CFdouble* rhoi, CFdouble* T) {}
  
  /**
   * Calculates the static pressure of the mixture
   * @param rho  density
   * @param temp temperature
   * @param tVec array of vibrational/electron temperatures
   */
  virtual CFdouble pressure(CFdouble& rho,
                            CFdouble& temp,
                            CFreal* tVec);


  /** Computes all transport coefficients for TCNEQ case
   *  @param rho  density
   *  @param temp temperature
   *  @param tVec array of vibrational/electron temperatures
   *  @param normConcGradients gradient of species mass fractions in one direction
   *  @param dynamic viscosity
   *  @param lambdaTrRo translational rotational thermal conductivity
   *  @param lambdaInt  internal (vibrational + electronic) thermal conductivity
   *  @param rhoUdiff   (1) Stefan-Maxwell model: rho*species diffusive velocity
   *                    (2) Fick, Fick+Ramshaw  : rowwise-ordered matrix of coefficients
   */
  virtual void transportCoeffNEQ(CFreal& temp,
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
   * @return total thermal conductivity (roto-translational, vibrational, electronic) for thermal equilibrium
   * This could be optional
   */
  virtual CFdouble lambdaNEQ(CFdouble& temp, CFdouble& pressure);

  /**
   * Calculates the thermal conductivity by conjugate gradient method method
   * given temperature and pressure
   * @param temp temperature
   * @param tVec vibrational/electronic temperatures
   * @param pressure pressure
   * @param lambdaTrRo roto-translational lambda (output)
   * @param lambdaInt  vibrational-electronic lambda (output)
   * @return total thermal conductivity (roto-translational, vibrational, electronic) for thermal NEQ
   */
  virtual void lambdaVibNEQ(CFreal& temp,
                            RealVector& tVec,
                            CFdouble& pressure,
                            CFreal& lambdaTrRo,
                            RealVector& lambdaInt);

  /**
   * Calculates the dynamic viscosity, given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   * @param tVec vibrational/electronic temperatures
   */
  virtual CFdouble eta(CFdouble& temp,
                       CFdouble& pressure,
                       CFreal* tVec);

  /**
   * Calculates the lambda viscosity, given temperature and pressure
   * @param temp temperature
   * @param pressure pressure
   */
  virtual CFdouble lambdaEQ(CFdouble& temp, CFdouble& pressure);

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
   * @param gamma specific heat ratio (output)
   * @param soundSpeed speed (output)
   * @param tVec vibrational/electronic temperatures (input)
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
  virtual void resetComposition(const RealVector& x);

  /**
   * Calculates the density, the enthalpy and the internal energy
   * This function is supplied for convenience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param tempV vibrational/electron temperature
   * @param pressure pressure
   * @param dhe array with density, enthalpy, energy. vibrational energies (one per vibrational equation) (output)
   * @param storeExtraData  flag telling if extra data have to be stored
   *                        (extra data are defined in src/Framework/PhysicalChemicalLibrary.hh)
   */
  virtual void setDensityEnthalpyEnergy(CFdouble& temp,
                                        RealVector& tVec,
                                        CFdouble& pressure,
                                        RealVector& dhe,
                                        bool storeExtraData);

  /**
   * Calculates the density, the enthalpy and the internal energy
   * This function is supplied for convenience, as it calls
   * density with electron temperature equal to temperature.
   * @param temp temperature
   * @param pressure pressure
   * @param dhe array with density, enthalpy, energy (output) for thermal equilibrium
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
   * Calculates the enthalpy in LTE conditions
   * at given temperature and pressure.
   * @param temp      temperature
   * @param pressure  pressure
   * @param intEnergy computed internal energy
   */
  virtual CFdouble enthalpy(CFdouble& temp,
                            CFdouble& pressure);

  /**
   * Calculates the speed of sound in
   * thermal equilibrium.
   * @param temp temperature
   * @param pressure pressure
   */
  virtual CFdouble soundSpeed(CFdouble& temp, CFdouble& pressure);

  /**
   * Sets the mole fractions of elements (nucleons) Xn for the ATDModel
   * environment (for variable elemental composition), this function
   * should be called before getting properties related to the elemental
   * fractions in case of LTE with Demixing
   * @param yn the RealVector of the mass fractions of elements
   */
  virtual void setElemFractions(const RealVector& yn);

  /**
   * Sets the mole fractions of elements (nucleons) Xn for the ATDModel
   * environment (for variable elemental composition) starting from the given
   * species mass fractions
   * @param yn the RealVector of the mass fractions of species
   */
  virtual void setElementXFromSpeciesY(const RealVector& ys);

  /**
   * Sets the species (molar) fractions. This function should be called before getting
   * thermodynamic quantities or transport properties.
   * @param ys the RealVector of the mass fractions of species
   */
  virtual void setSpeciesFractions(const RealVector& ys);

  /**
   * Sets the electron fractions in the mass composition according to charge neutrality.
   * @param ys the RealVector of the mass fractions of species
   */
  virtual void setElectronFraction(RealVector& ys);

  /**
   * Gets the species (molar) fractions.
   * @param ys the RealVector of the mass fractions of species (input)
   * @param xs the RealVector of the molar fractions of species (output)
   */
  virtual void getSpeciesMolarFractions(const RealVector& ys, RealVector& xs);

  /**
   * Gets the species (mass) fractions.
   * @param xs the RealVector of the molar fractions of species (input)
   * @param ys the RealVector of the mass fractions of species (output)
   */
  virtual void getSpeciesMassFractions(const RealVector& xs, RealVector& ys);

  /**
   * Gets the species mass fractions.
   * @param ys the RealVector of the mass fractions of species
   */
  virtual void getSpeciesMassFractions(RealVector& ys);

  /**
   * Returns the mass production/destruction terms [kg m^-3 s^-1] in CHEMICAL
   * NONEQUILIBRIUM based on Arrhenius's formula.
   * The analytical Jacobian matrix of the mass production terms can also be
   * computed for the variables rho_s, u, v,( w), T, T_v
   * @param pressure the mixture pressure
   * @param temp the mixture temperature
   * @param ys the species mass fractions
   * @param flg_jag the flag to compute the analytical Jacobian matrix
   * @param omega the mass production terms
   * @param jacobian the Jacobian matrix of the mass production terms (second step)
   */
  virtual void getMassProductionTerm(CFdouble& temp,
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
   * @param omegaRad ?
   */
  virtual void getSourceTermVT(CFdouble& temp,
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
  virtual void getSourceEE(CFdouble& temp,
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
   * @param omegaRad ?
   * @param jacobian the Jacobian matrix of the mass production terms
   */
  virtual void getSource(CFdouble& temp,
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
   * @param rhoUdiff   (1) Stefan-Maxwell model: rho*species diffusive velocity
   *                   (2) Fick, Fick+Ramshaw  : rowwise-ordered matrix of coefficients
   */
  virtual void getRhoUdiff(CFdouble& temp,
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
   * @param Dij diffusion coefficients
   * @param rhoUdiff rho*species diffusion velocities
   */
  virtual void getDij_fick(RealVector& dx,
                           CFdouble& pressure,
                           CFdouble& temperature,
                           RealMatrix& Dij,
                           RealVector& rhoUdiff);
  /**
   * Get the catalycity factor for N
   */
  virtual void getGammaN(CFreal& m_GN);

  /**
   * Get the catalycity factor for O
   */
  virtual void getGammaO(CFreal& m_GO);

  /**
   * Sets the species (molar) fractions. This function should be called before getting
   * thermodynamic quantities or transport properties.
   * @param xs the RealVector of the molar fractions of species
   */
  virtual void setSpeciesMolarFractions(const RealVector& xs);

  /**
   * Returns the total enthalpies per unit mass of species
   * @param tVec array of vibrational/electron temperatures
   * @param temp the mixture temperature
   * @param pressure the mixture pressure
   * @param hsTot species total enthalpy (output)
   * @param hsVib species vibrational enthalpy (output)
   * @param hsEl species electronic enthalpy (output)
   */
  virtual void getSpeciesTotEnthalpies(CFdouble& temp,
                                       RealVector& tVec,
                                       CFdouble& pressure,
                                       RealVector& hsTot,
                                       RealVector* hsVib,
                                       RealVector* hsEl);


private: // helper functions

  /**
   * This auxiliary class represents chemical reactions.
   *
   * @author Francisco Castro Rego
   *
   */
  class ChemReact{
  public:
    /**
     * Default Constructor
     */
    ChemReact();

    /**
     * Constructor for exchange reactions
     * @param stoichVec array with stoichiometric coefficients
     * @param chemCoefs array with coefficients to estimate reaction rates
     * @param chemID ID of chemistry mode
     * @param TempID Id of model for Tq
     */
    ChemReact(RealVector& stoichVec, RealVector& chemCoefs,CFuint& chemID,CFuint& TempID);

    /**
     * Constructor for dissociation reactions
     * @param stoichVec array with stoichiometric coefficients
     * @param chemCoefs array with coefficients to estimate reaction rates
     * @param dissPartners array indicating dissociation partners in reaction
     * @param chemID ID of chemistry mode
     * @param TempID Id of model for Tq
     */
    ChemReact(RealVector& stoichVec, RealVector& chemCoefs, std::vector<CFreal>& dissPartners,CFuint& chemID,CFuint& TempID);

    /**
     * Default Destructor
     */
    ~ChemReact();

    /**
     * Returns the reaction contribution to source term of species continuity equations.
     * @param temp the mixture temperature
     * @param tVec the vibrational temperature
     * @param pressure the mixture pressure
     * @param ys the species mass fractions
     * @param rho the mixture density
     * @param omega the mass production term
     */
    void omegaContribution(CFdouble& temperature, RealVector& tVec,CFdouble& pressure,CFdouble& rho,const RealVector& ys,const RealVector& mmasses,RealVector& omega);

    /**
     * Returns the forward reaction rate coefficient.
     * @param temperature the corrected temperature Tq
     */
    CFdouble Kf(CFdouble temperature);

    /**
     * Returns the backward reaction rate coefficient.
     * @param temperature the mixture temperature
     */
    CFdouble Kb(CFdouble temperature);

  private:

    ///  array with stoichiometric coefficients
    RealVector _stoichVec;

    ///  array with coefficients to estimate reaction rates
    RealVector _chemCoefs;

    ///  flag indicating dissociation reaction
    bool _flagDiss;

    /// ID of chemistry model: 0=park1 1=park2 2=park3 3=Dunn&Kang
    CFuint _chemID;

    /// Id of model for Tq: 0-Tq=sqrt(T*Tv) 1-Tq=T^0.7*Tv^0.3
    CFuint _TempID;

    ///  array indicating dissociation partners in reaction
    std::vector<CFreal> _dissPartners;
  };

  /**
   * Read mixture data
   */
  virtual void ReadDataMixture();

  /**
   * Read chemistry data
   */
  virtual void ReadDataChem();

  /**
   * Read species data
   * @param tVec array of vibrational/electron temperatures
   */
  virtual void ReadDataSpecies(CFint i);

  /**
   * Read thermo data
   * @param tVec array of vibrational/electron temperatures
   */
  virtual void ReadDataThermo(CFint i);

  /**
   * Set up the library
   */
  virtual void setLibrary();
  
  /**
   * Electronic-Vibrational energy of species over R (Candler)
   * @param tVib  vibrational/electron temperature
   * @param i index of species
   */
  CFdouble energyVibSpeciesOverR(CFdouble& tVib,CFint& i);

  /**
   * Total enthalpy on thermal equilibrium over R (Gnoffo)
   * @param temp temperature
   * @param i index of species
   */
  CFdouble enthalpyTeqOverRG(CFdouble& temp,CFint& i);

  /**
   * Total enthalpy on thermal equilibrium over R (McBride)
   * @param temp temperature
   * @param i index of species
   */
  CFdouble enthalpyTeqOverRMB(CFdouble& temp,CFint& i);

  /**
   * Electronic-Vibrational Cv of species over R (Candler)
   * @param tVib  vibrational/electron temperature
   * @param i index of species
   */
  CFdouble CvVibSpeciesOverR(CFdouble& tVib,CFint& i);

  /**
   * Cp on thermal equilibrium over R (Gnoffo)
   * @param temp temperature
   * @param i index of species
   */
  CFdouble CpTeqOverRG(CFdouble& temp,CFint& i);

  /**
   * Cp on thermal equilibrium over R (McBride)
   * @param temp temperature
   * @param i index of species
   */
  CFdouble CpTeqOverRMB(CFdouble& temp,CFint& i);

  /**
   * Electronic-Vibrational energy and Cv of species over R (Candler)
   * @param tVib  vibrational/electron temperature
   * @param i index of species
   * @param EvOR Electronic-Vibrational energy over R
   * @param CvvOR Electronic-Vibrational Cv over R
   */
  void energyVibSpeciesOverR(CFdouble& tVib,CFint& i,CFdouble& EvOR,CFdouble &CvvOR);

  /**
   * Total enthalpy and Cp on thermal equilibrium over R (Gnoffo)
   * @param temp temperature
   * @param i index of species
   * @param htotOR Total enthalpy over R
   * @param CpOR Specific heat in constant pressure over R
   */
  void enthalpyTeqOverRG(CFdouble& temp,CFint& i,CFdouble& htotOR,CFdouble &CpOR);

  /**
   * Total enthalpy and Cp on thermal equilibrium over R (McBride)
   * @param temp temperature
   * @param i index of species
   * @param htotOR Total enthalpy over R
   * @param CpOR Specific heat in constant pressure over R
   */
  void enthalpyTeqOverRMB(CFdouble& temp,CFint& i,CFdouble& htotOR,CFdouble &CpOR);

  /**
   * Read species involved in a chemical reaction from line
   * @param idx index of position in "line" separating reaction definition and coefficients
   * @param StchVec Vector containing Stoichiometric coefficients
   * @param getPartners flag indicating dissociation reaction
   * @param line line to be read
   */
  void readReactSpecies(CFint& idx,RealVector& StchVec,bool& getPartners,std::string& line);

  /**
   * Read coefficients to compute reaction rates from line
   * @param idx index of position in "line" separating reaction definition and coefficients
   * @param chemCoefs Coefficients to compute reaction rates
   * @param line line to be read
   */
  void readChemCoefs(CFint& idx,RealVector& chemCoefs,std::string& line);

  /**
   * Read partners of a dissociation reaction from line
   * @param dissPartners Vector with indexes of reaction partners
   * @param line line to be read
   */
  void readPartners(std::vector<CFreal>& dissPartners,std::string& line);

  /**
   * Returns the reaction contribution to source term of vibrational energy equation.
   * @param temp the mixture temperature
   * @param tVec the vibrational temperature
   * @param pressure the mixture pressure
   * @param ys the species mass fractions
   * @param rho the mixture density
   * @param omegav the vib energy source term
   */
  void ComputeOmegav(CFdouble& temperature, RealVector& tVec,CFdouble& pressure,CFdouble& rho,const RealVector& ys,RealVector& omega,RealVector& omegav);
    
protected: // data (first local arrays that are resized in setup(), then configurable options)

  ///  array with mass fractions
  RealVector _ys;

  ///  array with molar masses
  RealVector _mmasses;

  ///  array with formation enthalpies
  RealVector _hform;

  ///  array with characteristic temperature tauvs
  RealVector _tau_vs;

  ///  array Tets
  RealVector _Tets;

  ///  array Ss
  RealVector _Ss;

  ///  array As
  RealVector _As;

  ///  array Bs
  RealVector _Bs;

  ///  array Cs
  RealVector _Cs;

  ///  array with coefficients to compute enthalpy
  RealVector _thermoCoefs;

  ///  Moles of species over mass of mixture
  RealVector _massTtmS;

  ///  Auxiliary vector of stoichiometric coefficients
  RealVector _StchVecTemp;

  ///  Auxiliary vector of coefficients to compute reaction rates
  RealVector _chemCoefsTemp;

  ///  Species vibrational energy (temporary)
  RealVector _evibTemp;

  ///  Species vibrational energy evaluated at translational temperature (temporary)
  RealVector _evibAsterTemp;

  ///  Species vibrational heat capacity (temporary)
  RealVector _CvvsTemp;

  ///  Species vibrational heat capacity evaluated at translational temperature (temporary)
  RealVector _CvvsAsterTemp;

  ///  Mass production term (temporary)
  RealVector _omegaTemp;

  ///  Species viscosity (temporary)
  RealVector _NIUsTemp;

  ///  Matrix niusr to compute energy relaxation
  RealMatrix _NIUsr;

  ///  Matrix niusr to compute energy relaxation
  RealMatrix _Asr;
  
  /// mixture name
  std::string _mixtureName;

  /// thermo method name
  std::string _thermoName;

  /// chemistry model name
  std::string _chemName;

  /// array with names of species
  std::vector<std::string> _speciesNames;

  /// array with chemical reactions
  std::vector<ChemReact> _ChemReactArray;

  ///  example of array of IDs
  std::vector<CFuint> _exampleIDs;

  /// molecule IDs (example of array of IDs)
  std::vector<CFuint> _molecIDs;

  /// Constitutive atoms of species
  std::vector<CFuint> _CA;

  ///  Auxiliary vector of indexes of reaction partners
  std::vector<CFreal> _dissPartnersTemp;

  /// flag indicating molecules IDs
  std::vector<bool> _flagMoleculesIDs;

  /// Factor to reduce the stiffness of source terms <1
  CFdouble _factorOmega;

  /// Avogadro number
  CFdouble _NA;

  /// Boltzmann constant
  CFdouble _KB;

  /// pi
  CFdouble _PI;

  /// Lewis number
  CFdouble _Le;

  /// Moles over mass of mixture
  CFdouble _massTtmTt;

  /// Translational-Rotational specific heat for constant volume over R
  CFdouble _CvTrOR;

  /// Reference temperature
  CFdouble _Tref;

  /// Temperature at shock
  CFdouble _Tshk;

  /// Vibrational temperature at shock
  CFdouble _Tvsshk;

  /// ID of thermodynamic proprieties method: 0=Candler 1=Gnoffo 2=McBride
  CFuint _thermoID;

  /// ID of chemistry model: 0=park1 1=park2 2=park3 3=Dunn&Kang
  CFuint _chemID;

  /// dynamic viscosity model
  CFuint _etaModel;

  /// Id of model for Tq: 0-Tq=sqrt(T*Tv) 1-Tq=T^0.7*Tv^0.3
  CFuint _TempID;

  /// Id of model for Qt-v: 0-(Zhong's Report) 1-(Millican and White) 2-(Zhong's report without need for Tshk) 3-(Adapted from Mutation)
  CFuint _RelaxID;

  /// Id of model for Ss: 0-Species dependent 1-All 5000
  CFuint _TetsID;

  /// number of coefficients of the enthalpy model for each species
  CFint _numTC;

  /// number of coefficients of the enthalpy model for each species and range
  CFint _numTCR;

  /// number of chemical reactions
  CFint _numReac;

  }; // end of class ATDModelLibrary

//////////////////////////////////////////////////////////////////////////////

    } // namespace ATDModel

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ATDModel_ATDModelLibrary_hh
