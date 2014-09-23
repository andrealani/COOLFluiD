#ifndef COOLFluiD_Physics_ATDModel_CombustionModelLibrary_hh
#define COOLFluiD_Physics_ATDModel_CombustionModelLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "ATDModel/ATDModelLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ATDModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a CombustionModelLibrary.
 * For the order of the species in the mixture look *mix at ../Mutation2.0I/data/mixture/
 *
 * @author Alessandro Mazzetti
 *
 */
class CombustionModelLibrary : public ATDModelLibrary {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  CombustionModelLibrary(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CombustionModelLibrary();

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Unsetups the data of the library
   */
  virtual void unsetup();
  
  /**
   * Set the IDs of the molecules in the mixture
   */
  virtual void setMoleculesIDs(std::vector<CFuint>& v);
  
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
                                 CFreal& eta,
                                 CFreal& lambdaTrRo,
                                 RealVector& lambdaInt,
                                 RealVector& rhoUdiff);

  /**
   * Calculates the thermal conductivity given temperature, from Svehla 1995.
   * @param temp temperature
   * @param pressure pressure
   * @return total thermal conductivity (roto-translational, vibrational, electronic) for thermal equilibrium
   * This could be optional
   */
  virtual CFdouble lambdaNEQ(CFdouble& temp, CFdouble& pressure);//,CFreal* tVec);
  //Before it was (CFdouble& temp, CFdouble& pressure);

  /**
   * Calculates the dynamic viscosity, given temperature, from Svehla 1995.
   * @param temp temperature
   * @param pressure pressure
   * @param tVec vibrational/electronic temperatures
   */
  virtual CFdouble eta(CFdouble& temp, 
		       CFdouble& pressure, 
		       CFreal* tVec);
  
   /**
   * Species Cp calculation with polynomial interpolation
   * @param temp temperature
   * @param i index of species
   */
  virtual CFdouble CpPoly(CFdouble& temp);//,CFint& i);
  
   /**
   * Species enthalpy calculation with polynomial interpolation
   * @param temp temperature
   * @param i index of species
   */
  virtual void EnthPoly(CFdouble& temp,CFint& i);//,CFdouble& _hsP);

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
                                        CFdouble& soundSpeed,//,
                                        RealVector* tVec);

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
                                        //RealVector& tVec,
                                        CFdouble& pressure,
                                        RealVector& dhe);//,
                                        //bool storeExtraData);
  
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
                           CFreal& Diff_coeff,
                           RealMatrix& Dij,
                           RealVector& rhoUdiff);
  
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
/**
   * Returns the omega mass production terms
   * @param temp the mixture temperature
   * @param rho the mixture density
   * @param ys the species mass fractions
   * @param mmasses the species mole fraction
   * @param omega the mass production term vector, ordered by species
   */
  virtual void omegaContribution(CFdouble& temperature,
				 RealVector& tVec,
				 CFdouble& pressure,
				 CFdouble& rho,
				 const RealVector& ys,
				 const RealVector& mmasses,
				 RealVector& omega);
                                 

private: // helper functions

  /// if you need to use ChemReact you have to use the full name ATDModelLibrary::ChemReact
  
  /**
   * Read mixture data
   */
  virtual void ReadDataMixture();

  /**
   * Read chemistry data
   */
//  virtual void ReadDataChem();
  
    /**
   * Read species data
   * @param tVec array of vibrational/electron temperatures
   */
  virtual void ReadDataSpecies(CFint i);

  /**
   * Read thermo data
   * @param tVec array of vibrational/electron temperatures
   */
  virtual void ReadDataThermo(CFint i);// {/*check if it is not needed at all*/ }
  
  /**
   * Set up the library
   */
  virtual void setLibrary();
  
protected: // data (first local arrays that are resized in setup(), then configurable options)


// _AsMu,_AsK, _BsMu, _BsK and so on are vectors of length "Number of species - 1" (butadiene is calculated otherwise and
// it is the first member of species viscosity array); everyone contains one out of four coefficients for
// species viscosity calculation from Svehla et al., 1995, in the LOWER temperature range.
// Higher temperature range coefficients (above 1000K, usually) are calculated run-time with a numerical correction value
// calculated as ratio e.g. of _AsMu_HigherT_Range/_AsMu_LowerT_Range

///  array of molar fractions
RealVector _xs;

///  array _AsMu1;
RealVector _AsMu1;

///  array _BsMu1;
RealVector _BsMu1;

///  array _CsMu1;
RealVector _CsMu1;

///  array _DsMu1;
RealVector _DsMu1;

///  array _AsMu2;
RealVector _AsMu2;

///  array _BsMu2;
RealVector _BsMu2;

///  array _CsMu2;
RealVector _CsMu2;

///  array _DsMu2;
RealVector _DsMu2;

///  array _AsK1;
RealVector _AsK1;

///  array _BsK1;
RealVector _BsK1;

///  array _CsK1;
RealVector _CsK1;

///  array _DsK1;
RealVector _DsK1;

///  array _AsK2;
RealVector _AsK2;

///  array _BsK2;
RealVector _BsK2;

///  array _CsK2;
RealVector _CsK2;

///  array _DsK2;
RealVector _DsK2;

//_a(i)cp(j) and _a(i)h(j) are the coefficients of polynomial interpolation for
// specific heat and enthalpy calculation; each array contains coefficients accompanying
// highest to lowest power elevation of temperature.
// (i) is 1 or 2, expressing 200-1000 or 1000-6000 K temperature range
// (j) is 1 to 5 for cp, and 1 to 6 for enthalpy
// each vector contains that "order" coefficients for all species in the order 
// expressed in _speciesNames[i] vector.

///  array _a1cp1;
RealVector _a1cp1;

///  array _a2cp1;
RealVector _a2cp1;

///  array _a3cp1;
RealVector _a3cp1;

///  array _a4cp1;
RealVector _a4cp1;

///  array _a5cp1;
RealVector _a5cp1;

///  array _a1cp2;
RealVector _a1cp2;

///  array _a2cp2;
RealVector _a2cp2;

///  array _a3cp2;
RealVector _a3cp2;

///  array _a4cp2;
RealVector _a4cp2;

///  array _a5cp2;
RealVector _a5cp2;

///  array _a1h1;
RealVector _a1h1;

///  array _a2h1;
RealVector _a2h1;

///  array _a3h1;
RealVector _a3h1;

///  array _a4h1;
RealVector _a4h1;

///  array _a5h1;
RealVector _a5h1;

///  array _a6h1;
RealVector _a6h1;

///  array _a1h2;
RealVector _a1h2;

///  array _a2h2;
RealVector _a2h2;

///  array _a3h2;
RealVector _a3h2;

///  array _a4h2;
RealVector _a4h2;

///  array _a5h2;
RealVector _a5h2;

///  array _a6h2;
RealVector _a6h2;

///  Species Thermal Conductivity (temporary)
RealVector _KappasTemp;

///  Species Specific Heat (temporary)
RealVector _CpTemp;

///  Species Specific Heat (temporary)
RealVector _hsTemp;

///  Species Specific Heat (temporary)
RealVector _hsP;

/// Mass prouction term (temporary)
RealVector _omegaTempDot;

/// Species total mass (molar mass over mass fraction ratio)
//RealVector _MtotS;

/// Acentric Factor
CFdouble _AF;

/// Critical Temperature
CFdouble _Tc;

/// Adimensional Temperature
CFdouble _Tr;

/// Critical Volume
CFdouble _Vc;

/// Epsilon over kb
CFdouble _EpsKb;

/// Correction Factor
CFdouble _Fc;

/// Collision Factor Omega_v
CFdouble _Ov;

/// Adimensional temperature
CFdouble _Tstar;

/// mean molar mass
CFdouble _meanmass;

// /// Specific Heat, calculated with polynomial interpolation
//CFdouble _CpP;

/// Specific Heat for C4H6
CFdouble _CpC4H6;

/// Universal gas constant for cal/(mol K) use in chemical reaction rate
CFdouble _Rgascal;

/// Forward reaction coefficient 1, Venkateswaran chem. model
CFdouble _kfw1;

/// Forward reaction coefficient 2, Venkateswaran chem. model
CFdouble _kfw2;

/// Backward reaction coefficient 2, Venkateswaran chem. model
CFdouble _kbw2;

/// Reaction coefficient 1, Jones-Lindstedt chem. model
CFdouble _ka;

/// Reaction coefficient 2, Jones-Lindstedt chem. model
CFdouble _kb;

/// Reaction coefficient 3, Jones-Lindstedt chem. model
CFdouble _kc;

/// Reaction coefficient 4, Jones-Lindstedt chem. model
CFdouble _kd;

/// Reaction coefficient 5, Jones-Lindstedt chem. model
CFdouble _ke;

/// Reaction coefficient 6, Jones-Lindstedt chem. model
CFdouble _kf;

/// Reaction rate 1, Jones-Lindstedt chem. model
CFdouble _Ra;

/// Reaction rate 2, Jones-Lindstedt chem. model
CFdouble _Rb;

/// Reaction rate 3, Jones-Lindstedt chem. model
CFdouble _Rc;

/// Reaction rate 4, Jones-Lindstedt chem. model
CFdouble _Rd;

/// Reaction rate 5, Jones-Lindstedt chem. model
CFdouble _Re;

/// Reaction rate 6, Jones-Lindstedt chem. model
CFdouble _Rf;

/// Mass Diffusion coefficient, constant for all species
CFdouble _D;

/// Mixture total mass (sum for all species of molar mass over mass fraction ratio)
//CFdouble _Mtot;

}; // end of class CombustionModelLibrary

//////////////////////////////////////////////////////////////////////////////

    } // namespace ATDModel

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ATDModel_CombustionModelLibrary_hh
