#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCat_nad_hh 
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCat_nad_hh 

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      template <class BASEVS> class MultiScalarVarSet;
      class Euler2DVarSet;
    }
  }

  
  namespace Framework {
    class PhysicalChemicalLibrary;
    class CatalycityModel;
  }
  
  namespace MathTools {
    class MatrixInverter;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
 /**
  * This class implements a catalytic boundary condition
  * This is still a very draft version, that is using
  * A diffusive flux of Fick (with the Di constant)
  * do not consider ions, only one temperature. It is working only for air5
  * @author Nadege Villedieu
  *
  */
   
class NoSlipWallIsothermalNSrvtCat_nad : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSrvtCat_nad(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSrvtCat_nad();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

 private: // data
  /// Function that computes the mass fraction of each spicies
  /// @param innerState is the state of the inner cell
  /// @param ghostState is the state of the ghost cell
  /// It computes the massic fraction and the total densities

  void compute_mass_fraction(Framework::State& innerState, Framework::State& ghostState, RealVector& m_yb, RealVector& m_yi, CFreal& m_rhob, CFreal& m_rhoi);

  /// Function that do the transformation from massique
  /// to molar fraction
  void mass_frac_to_mol_frac(RealVector& m_yb, RealVector& m_yi, RealVector& m_xb, RealVector& m_xi);

  
  /** This function is computing the rate production at the wall
   * Returns the total enthalpies per unit mass of species
   * @param m_rhob the density at the wall
   * @param m_yb  the massic fractions
   * @param m_Tw  temperature of the wall
   * @param m_gamma catalicity factor    
   * TODO Now, gamma is the same for all equations and spieces
   * This needs to be improved
    * @param m_nu matrix of the destructions
    * @param m_muN matrix of the production
    * TODO Now, this is working for air5 considering that at
    * the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO  **/
  void getWallRateProd(CFreal& m_rhob,
                       RealVector& m_yb,
                       CFreal& m_Tw,
                       RealMatrix& m_nu,
                       RealMatrix& m_muN,
                       RealMatrix& m_muO,
                       RealMatrix& m_muN2,
                       RealMatrix& m_muNO,
                       RealMatrix& m_muO2,
                       RealVector& m_omegawall);
  
  /// This function transform the molar fraction in partial density
  ///@param m_xb molar fraction of the ghost cell
  ///@param m_p pressure of the inner cell
  ///@param m_Tw Temperature at the wall
  void mol_frac_to_part_dens(RealVector& m_xb, CFreal& m_p, CFreal& m_Tw, RealVector& m_partrho);
  void mol_frac_to_mass_frac(RealVector& m_xb, RealVector& m_yb);
  void getMolarMass(RealVector &m_xp,CFreal &m_mmt);
  /// The socket to use in this strategy for the distance to the wall
  Framework::DataSocketSink<CFreal> socket_walldistance;


  /// temporary storage of the wall distance
  CFreal m_walldistance ;

   /// Temperature at the wall
   CFreal m_Tw;

  /// Number of spices
  CFuint m_nbSpecies;

  /// Massique fractions at the interior
  RealVector m_yi;

  /// Molar fractions at the interior
  RealVector m_xi;

  /// Masique fractions at the wall
  RealVector m_yb;

  /// perturbed Masique fractions at the wall
  RealVector m_yp;

  /// Molar fractions at the wall
  RealVector m_xb;

  /// Molar fractions at the wall that we perturbe to compute 
  /// the diffusive flux
  RealVector m_xp;

  /// Vector of the molar massestha we perturbe to compute 
  /// the diffusive flux
  RealVector m_dxp;

  /// Vector of the molar masses
  RealVector m_mm;

  /// Vector of the derivative of the xb
  RealVector m_dx;

  /// Matrice of the diffusion coefficients (since we consider them vconstant here, it is not used)
  RealMatrix m_Dbin;

  /// Vector of the diffusive flux
  RealVector m_Diff_flux;
  
  /// Vector of the diffusive flux when we perturbe
  RealVector m_Diff_fluxp;

  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DVarSet> > _varSet;

  /// Total density at  the wall (in the ghost cell)
  CFreal m_rhob;

  /// Total density at  the wall
  CFreal m_rhoi;

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// Number of reactions catalysed by the wall
  CFreal m_nr;

  /// Gamma of the wall 
  /// TODO for the moment it is using the same gamma for all spicies
  /// And all reactions
  CFreal m_GammaN;
  CFreal m_GammaO;

   /// vector conatining the nu of each reaction
   /// TODO for the moment we suppose that we use air5 and
   /// that the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO
   RealMatrix m_nu;
  
  /// vector conatining the mu for N of each reactions
  /// TODO for the moment we suppose that we use air5 and
  /// that the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO
  RealMatrix m_muN;

  /// vector conatining the mu for O of each reactions
  /// TODO for the moment we suppose that we use air5 and
  /// that the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO
  RealMatrix m_muO;

  /// vector conatining the mu for N2 of each reactions
  /// TODO for the moment we suppose that we use air5 and
  /// that the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO
  RealMatrix m_muN2;

  /// vector conatining the mu for NO of each reactions
  /// TODO for the moment we suppose that we use air5 and
  /// that the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO
  RealMatrix m_muNO;

  /// vector conatining the mu for O2 of each reactions
  /// TODO for the moment we suppose that we use air5 and
  /// that the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO
  RealMatrix m_muO2;
  
  /// Prodiction rate at the wall
  RealVector m_omega;
  RealVector m_omegap;

  /// Vector of the impinging fluxes
  RealVector m_mcal;

  /// Number of Newton loop
  CFuint m_nl;

  /// Second hand right-hand side 
  RealVector m_b;
  
/// Second hand right-hand side 
  RealVector m_zero;
   /// Matrix of the derivative of the Diffusion flux
   /// minus the wall production wall
   RealMatrix m_a;

  /// Inverse matrix of the derivative of the Diffusion flux
   /// minus the wall production wall
   RealMatrix m_inva;

  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter> m_inverter;

  /// ID of the temperature
  CFuint m_tempID;

  /// Partial densities
  RealVector m_partrho;
  
  /// temporary ghost rhoi
  RealVector m_rhoG;

}; // end of class NoSlipWallIsothermalNSrvtCat-nad

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCat_nad_hh
