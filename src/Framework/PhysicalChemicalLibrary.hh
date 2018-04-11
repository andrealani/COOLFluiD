// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PhysicalChemicalLibrary_hh
#define COOLFluiD_Framework_PhysicalChemicalLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalPropertyLibrary.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a physico-chemical library
/// @author Andrea Lani
class Framework_API PhysicalChemicalLibrary : public Framework::PhysicalPropertyLibrary {
public:

  class Framework_API ExtraData {
  public:

    ExtraData(){needsBarVars = false;}

    ~ExtraData(){}

    /// species total entalpy
    RealVector enthalpyTt;

    /// species translational-rotational energies
    RealVector energyTr;

    /// species vibrational energies
    RealVector energyVib; 

    /// species vibrational specific heat
    RealVector cpVib;  

    /// species electronic specific heat
    //  for free electrons the translational specific heat is stored
    RealVector cpElec;

    /// species electronic energy
    //  for free electrons the translational energy is stored
    RealVector eElec;

    /// species formation enthalpies
    RealVector enthalpyForm;

    /// derivative of density*energy with partial densities
    RealVector dRhoEdRhoi;

    /// derivative of density*energy with partial densities
    RealVector dRhoEvdRhoi;

    /// derivative of internal energy with T
    CFreal dEdT;

    /// derivative of internal energy with T
    CFreal dHdT;

    /// derivative of vib energy with Tv
    CFreal dEvTv;

    /// molecules in the vib energy equations
    std::vector<CFuint> list;

    /// Bar variables for Roe-Vinokur linearization
    bool needsBarVars;

    /// Vector storing dpdRhoi and dpdRhoE at bar conditions (Roe-Vinokur linearization)
    /// Components (0,Ns-1) contain dpdRhoi, while component Ns contains kappa = dpdRhoE
    RealVector dP_Bar;

  };

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor without arguments
  PhysicalChemicalLibrary(const std::string& name);

  /// Default destructor
  virtual ~PhysicalChemicalLibrary();

  /// Configures this configurable object.
  virtual void configure ( Config::ConfigArgs& args );

  /// Set up the number of vibrational temperatures
  virtual void setNbTempVib(CFuint nbTv)
  {
    _nbTvib = nbTv;
  }

  /// Set the number of electron temperature (0 or 1)
  virtual void setNbTe(CFuint nbTe)
  {
    _nbTe = nbTe;
  }
  
  /// Get the number of vibrational temperatures
  CFuint getNbTe() const
  {
    return _nbTe;
  }

  /// Get the ID of the variable entry that gets the electronic energy
  CFint getElectrEnergyID() const
  {
    return _electrEnergyID;
  }

  /// Get the number of vibrational temperatures
  CFuint getNbTempVib() const
  {
    return _nbTvib;
  }

  /// Get the number of species
  CFint getNbSpecies() const
  {
    return _NS;
  }

  /// Get the number of elements
  CFint getNbElements()
  {
    return _NC;
  }

  /// Get the constant of perfect gases in J/(mol*K)
  CFdouble getRgas() const
  {
    return _Rgas;
  }

  /// Get the IDs corresponding to the molecules for each vibrational energy equation
  const std::vector<CFuint>& getMolecule2EqIDs() const
  {
    cf_assert(_molecule2EqIDs.size() > 0);
    return _molecule2EqIDs;
  }

  /// Get the electron pressure
  CFdouble getElectronPress() const
  {
    return _electronPress;
  }

  /// Get the atomicity coeficient
  Common::SafePtr<RealVector> getAtomicityCoeff() 
  {
    return &_atomicityCoeff;
  }

  /// Compute and get the electron pressure
  virtual CFdouble electronPressure(CFreal rhoE,
  			    CFreal tempE) = 0;

  /// Get the molar masses
  virtual void getMolarMasses(RealVector& mm) = 0;

  /// Get Catalycity Factor for N
  virtual void getGammaN(CFreal& m_GN) 
  {
    throw Common::NotImplementedException(FromHere(),"PhysicalChemicalLibrary::getGammaN()");
  }
  
  /// Get Catalycity Factor for O
  virtual void getGammaO(CFreal& m_GO) 
  {
    throw Common::NotImplementedException(FromHere(),"PhysicalChemicalLibrary::getGammaO()");
  }
  
  /// Get the translational-rotational cv
  /// @pre it assumes that the mass fractions
  ///      have been already set
  virtual CFdouble getCvTr() const = 0;

  /// Get the translational-rotational cv
  /// @pre it assumes that the mass fractions
  ///      have been already set
  virtual CFdouble getMMass() const = 0;

  /// Set the constant of gases in J/(Kg*K)
  virtual void setRiGas(RealVector& Ri) = 0;

  /// Set the IDs of the molecules in the mixture
  virtual void setMoleculesIDs(std::vector<CFuint>& v) = 0;
               
  /// return true if there is some electron in the mixture (eg: true for air11 and false for air5)
  bool presenceElectron() const
  {
    return _hasElectrons;
  }
  
  /// Set the thermodynamic state (temperature and pressure)
  /// @param species partial densities
  /// @param mixture temperature
  virtual void setState(CFdouble* rhoi, CFdouble* T) = 0;
  
  /// Calculates the static pressure of the mixture
  /// @param rho  density
  /// @param temp temperature
  virtual CFdouble pressure(CFdouble& rho,
			    CFdouble& temp,
			    CFreal* tVec) = 0;
  
  /// Computes all transport coefficients for TCNEQ case
  virtual void transportCoeffNEQ(CFreal& temp, 
				 CFdouble& pressure,
				 CFreal* tVec, 
				 RealVector& normConcGradients,
				 RealVector& normTempGradients,
				 CFreal& eta,
				 CFreal& lambdaTrRo, 
				 RealVector& lambdaInt,
				 RealVector& rhoUdiff) = 0;
  
  /// Calculates the thermal conductivity by conjugate gradient method method
  /// given temperature and pressure
  /// @param temp temperature
  /// @param pressure pressure
  virtual CFdouble lambdaNEQ(CFdouble& temp, CFdouble& pressure) = 0;

  /// Calculates the thermal conductivity by conjugate gradient method method
  /// given temperature and pressure
  /// @param temp temperature
  /// @param pressure pressure
  virtual void lambdaVibNEQ(CFreal& temp,
      RealVector& tVec,
      CFdouble& pressure,
      CFreal& lambdaTrRo,
      RealVector& lambdaVib)  = 0;

  /// Calculates the dynamic viscosity, given temperature and pressure
  /// @param temp temperature
  /// @param pressure pressure
  virtual CFdouble eta(CFdouble& temp,
  	       CFdouble& pressure,
  	       CFreal* tVec) = 0;

  /// Calculates the lambda viscosity, given temperature and pressure
  /// @param temp temperature
  /// @param pressure pressure
  virtual CFdouble lambdaEQ(CFdouble& temp,
      CFdouble& pressure) = 0;

  /// Calculates the electrical conductivity given temperature and pressure
  /// @param temp temperature
  /// @param pressure pressure
  virtual CFdouble sigma(CFdouble& temp,
  		 CFdouble& pressure,
  		 CFreal* tVec) = 0;

  /// Calculates the specific heat ratio and the speed of sound in
  /// thermal equilibrium.
  /// @param temp temperature
  /// @param pressure pressure
  /// @param rho density
  /// @param gamma specific heat ratio
  /// @param soundSpeed speed
  virtual void gammaAndSoundSpeed(CFdouble& temp,
  			  CFdouble& pressure,
  			  CFdouble& rho,
  			  CFdouble& gamma,
  			  CFdouble& soundSpeed) = 0;

  /// Calculates the frozen speed of sound in
  /// thermal equilibrium
  /// @param temp temperature
  /// @param pressure pressure
  /// @param rho density
  virtual void frozenGammaAndSoundSpeed(CFdouble& temp,
    CFdouble& pressure,
    CFdouble& rho,
    CFdouble& gamma,
    CFdouble& soundSpeed,
    RealVector* tVec) = 0;

  /// Calculates the composition given temperature and pressure.
  /// @param temp temperature
  /// @param pressure pressure
  /// @param x composition array (one component for each species)
  /// @pre this function ALWAYS computes the composition
  ///      (no lookup table based results...)
  virtual void setComposition(CFdouble& temp,
  		      CFdouble& pressure,
  		      RealVector* x = CFNULL) = 0;

  /// Reset the composition
  virtual void resetComposition(const RealVector& x) = 0;

  /// Calculates the density given temperature and pressure.
  /// This function is supplied for convinience, as it calls
  /// density with electron temperature equal to temperature.
  /// @param temp temperature
  /// @param pressure pressure
  virtual CFdouble density(CFdouble& temp,
  		   CFdouble& pressure,
  		   CFreal* tVec = CFNULL) = 0;

  /// Calculates the density, the enthalpy and the internal energy
  /// This function is supplied for convenience, as it calls
  /// density with electron temperature equal to temperature.
  /// @param temp temperature
  /// @param tempV vibrational temperature
  /// @param pressure pressure
  /// @param dhe array with density, enthalpy, energy. vibrational energies
  /// @param storeExtraData  flag telling if extra data have to be stored
  virtual void setDensityEnthalpyEnergy(CFdouble& temp,
    RealVector& tVec,
    CFdouble& pressure,
    RealVector& dhe,
    bool storeExtraData = false) = 0;

  /// Calculates the density, the enthalpy and the internal energy
  /// This function is supplied for convenience, as it calls
  /// density with electron temperature equal to temperature.
  /// @param temp temperature
  /// @param pressure pressure
  virtual void setDensityEnthalpyEnergy(CFdouble& temp,
    CFdouble& pressure,
    RealVector& dhe) = 0;

  /// Get some extra data like the species translational-rotational
  /// energies or some partial derivatives that have been previously
  /// computed and stored
  Common::SafePtr<ExtraData> getExtraData()
  {
    return &_extraData;
  } 

  /// Sets the mole fractions of elements (nucleons) Xn for the Mutation
  /// environment (for variable elemental composition), this function
  /// should be called before getting properties related to the elemental
  /// fractions in case of LTE with Demixing
  /// @param yn the RealVector of the mass fractions of elements
  virtual void setElemFractions(const RealVector& yn) = 0;

  /// Sets the mole fractions of elements (nucleons) Xn for the Mutation
  /// environment (for variable elemental composition) starting from the given
  /// species mass fractions
  /// @param yn the RealVector of the mass fractions of species
  virtual void setElementXFromSpeciesY(const RealVector& ys) = 0;

  /// Sets the species (molar) fractions. This function should be called before getting
  /// thermodynamic quantities or transport properties.
  /// @param ys the RealVector of the mass fractions of species
  virtual void setSpeciesFractions(const RealVector& ys) = 0;

  /// Sets the species (molar) fractions. This function should be called before getting
  /// thermodynamic quantities or transport properties.
  /// @param ys the RealVector of the mass fractions of species
  virtual void setSpeciesMolarFractions(const RealVector& xs) = 0;
  
  /// Sets the electron fractions in the mass composition according to charge
  /// neutrality.
  /// @param ys the RealVector of the mass fractions of species
  virtual void setElectronFraction(RealVector& ys) = 0;

  /// Gets the species (molar) fractions.
  /// @param ys the RealVector of the mass fractions of species
  /// @param xs the RealVector of the molar fractions of species
  virtual void getSpeciesMolarFractions(const RealVector& ys, RealVector& xs) = 0;

  /// Gets the species (mass) fractions.
  /// @param xs the RealVector of the molar fractions of species
  /// @param ys the RealVector of the mass fractions of species
  virtual void getSpeciesMassFractions(const RealVector& xs, RealVector& ys) = 0;

  /// Gets the species mass fractions.
  /// @param ys the RealVector of the mass fractions of species
  virtual void getSpeciesMassFractions(RealVector& ys) = 0;

  /// Returns the transport coefficients for the diffusive fluxes used in
  /// computations with LTE and Demixing
  /// @param temp the mixture temperature
  /// @param pressure the mixture pressure
  /// @param lambda the Butler and Brokaw thermal reactive conductivity
  /// @param lambdacor the demixing thermal reactive conductivity
  /// @param lambdael the elmental heat transfer coefficients
  /// @param eldifcoef elemental multicomponent difussion coefficients
  ///        times mixture density
  /// @param eltdifcoef the elemental thermal demixing coefficients
  ///        times mixture density
  virtual void getTransportCoefs(CFdouble& temp,
				 CFdouble& pressure,
				 CFdouble& lambda,
				 CFdouble& lambdacor,
				 RealVector& lambdael,
				 RealMatrix& eldifcoef,
				 RealVector& eltdifcoef)
  {
    throw Common::NotImplementedException(FromHere(),"PhysicalChemicalLibrary::getTransportCoefs()");
  }
  
  /// Returns the mass production/destruction terms [kg m^-3 s^-1] in CHEMICAL
  /// NONEQUILIBRIUM based on Arrhenius's formula.
  /// The analytical Jacobian matrix of the mass production terms can also be
  /// computed for the variables p, u, v,( w,) ys
  /// @param pressure the mixture pressure
  /// @param temp the mixture temperature
  /// @param ys the species mass fractions
  /// @param flg_jag the flag to compute the analytical Jacobian matrix
  /// @param omega the mass production terms
  /// @param jacobian the Jacobian matrix of the mass production terms
  virtual void getMassProductionTerm(CFdouble& temp,
  			     RealVector& tVec,
  			     CFdouble& pressure,
  			     CFdouble& rho,
  			     const RealVector& ys,
  			     bool flagJac,
  			     RealVector& omega,
  			     RealMatrix& jacobian) = 0;

  /// Returns the source term for the vibrational relaxation with VT transfer
  /// @param temp the mixture temperature
  /// @param tVec the vibrational temperature
  /// @param pressure the mixture pressure
  /// @param rho the mixture density
  /// @param omegav the source term
  virtual void getSourceTermVT(CFdouble& temp,
  		       RealVector& tVec,
  		       CFdouble& pressure,
  		       CFdouble& rho,
  		       RealVector& omegav,
                         CFdouble& omegaRad) = 0;

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
   virtual void getSource(CFdouble& temp,
			  RealVector& tVec,
			  CFdouble& pressure,
			  CFdouble& rho,
			  const RealVector& ys,
			  bool flagJac,
			  RealVector& omega,
			  RealVector& omegav,
			  CFdouble& omegaRad,
			  RealMatrix& jacobian) = 0;
  
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
			   CFdouble& omegaEE) = 0;
  
 
  /// Returns the diffusion velocities of species multiplied by the species
  /// densities for nonequilibrium computations
  /// @param temp the mixture temperature
  /// @param pressure the mixture pressure
  /// @param normConcGradients the cell normal gradients of species mass fractions
  virtual void getRhoUdiff(CFdouble& temp,
			   CFdouble& pressure,
			   RealVector& normConcGradients,
			   RealVector& normTempGradients,
			   CFreal* tVec,
			   RealVector& rhoUdiff,
			   bool fast = false) = 0;

  /// Returns the diffusion velocities of species multiplied by the species
  /// densities for nonequilibrium computations
  /// @param temp the mixture temperature
  /// @param pressure the mixture pressure
  
  virtual void getDij_fick(RealVector& dx,
			   CFdouble& pressure,
			   CFdouble& temperature,
			   RealMatrix& Dij,
			   RealVector& rhoUdiff)
  {
      throw Common::NotImplementedException(
          FromHere(),"PhysicalChemicalLibrary::getDij_fick()");
  }
  
  /// Returns the total enthalpies per unit mass of species
  /// @param temp the mixture temperature
  /// @param pressure the mixture pressure 
  virtual void getSpeciesTotEnthalpies(CFdouble& temp,
				       RealVector& tVec,
				       CFdouble& pressure,
				       RealVector& hsTot,
				       RealVector* hsVib = CFNULL,
				       RealVector* hsEl = CFNULL) = 0;
  
  /// Temperature of free electrons
  CFdouble getTe(CFdouble temp, CFreal* tVec)
  {
    // tVec can be CFNULL
    if (tVec == CFNULL) return temp;
    if (_electrEnergyID < 0) return temp;
    cf_assert(_electrEnergyID < static_cast<CFint>((_nbTvib + getNbTe())));
    return (presenceElectron()) ? std::min(tVec[_electrEnergyID],(CFreal)_maxTe) : temp;
  }
  
protected:
  
  /// number of (types of) species
  int _NS;

  /// Number of (types of) elements
  int _NC;

  /// number of vibrational temperatures
  int _nbTvib;

  /// number of free electron temperatures (0 or 1)
  int _nbTe;
  
  /// flag telling if the mixture is ionized
  bool _hasElectrons;
  
  /// constant of perfect gases
  CFdouble _Rgas;
  
  /// electron pressure
  CFdouble _electronPress;
 
  /// mixture formation enthalphy at T=0K
  CFdouble m_H0;
  
  /// species formation enthalphy at T=0K
  RealVector m_vecH0;
  
  // some additional data
  ExtraData _extraData;
  
  /// atomicity coefficient
  RealVector _atomicityCoeff;
 
  /// IDs corresponding to the molecules for each vibrational energy equation
  std::vector<CFuint> _molecule2EqIDs;

  /// ID of the variable entry that gets the electronic energy
  CFint _electrEnergyID;

  /// freeze the chemistry
  bool  _freezeChemistry;
  
  /// shift the formation enthalpy to have H(T=0K)=0
  bool m_shiftHO;
  
  /// Max value for Te
  CFdouble _maxTe;
  
}; // end of class PhysicalChemicalLibrary

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PhysicalChemicalLibrary_hh
