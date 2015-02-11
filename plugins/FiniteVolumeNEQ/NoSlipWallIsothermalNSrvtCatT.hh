#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCatT_hh 
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCatT_hh 

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "MathTools/FunctionParser.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
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
  *
  * @author Nadege Villedieu
  * @author Andrea Lani
  *
  */
template <class MODEL>
class NoSlipWallIsothermalNSrvtCatT : public NoSlipWallIsothermalNSrvt<MODEL> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSrvtCatT(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSrvtCatT();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /// Configures the command.
  virtual void configure ( Config::ConfigArgs& args );

protected:
  
  /**
   * Apply boundary condition on the given ghost state
   * @param innerState state vector associated to the internal cell
   * @param ghostState ghost state vector associated to the ghost cell
   * @param 
   */
  virtual void setGhostStateImpl(const Framework::State& innerState, 
				 Framework::State& ghostState);
  
private: // data
  /// Function that computes the mass fraction of each species
  /// @param innerState is the state of the inner cell
  /// @param ghostState is the state of the ghost cell
  /// It computes the massic fraction and the total densities

  void compute_mass_fraction(const Framework::State& innerState, 
			     const Framework::State& ghostState, 
			     RealVector& yb, RealVector& yi, 
			     CFreal& rhob, CFreal& rhoi);

  /** This function is computing the rate production at the wall
   * Returns the total enthalpies per unit mass of species
   * @param press pressure at the wall
   * @param yb    the mass fractions
   * @param Tw    temperature of the wall
   * @param gamma catalicity factor    
   * TODO Now, gamma is the same for all equations and spieces
   * This needs to be improved
   * @param nu matrix of the destructions
   * @param muN matrix of the production
   * TODO Now, this is working for air5 considering that at
   * the catalised reactions are N+N -> N2; O+O -> O2; N+O -> NO  **/
  void getWallRateProd(const CFreal press,
                       const RealVector& yb,
                       const CFreal Tw,
                       const RealMatrix& nu,
		       const RealMatrix& muN,
                       const RealMatrix& muO,
                       const RealMatrix& muN2,
                       const RealMatrix& muNO,
                       const RealMatrix& muO2,
		       RealVector& omegawall);
  
  /// This function transform the molar fraction in partial density
  ///@param xb molar fraction of the ghost cell
  ///@param p pressure of the inner cell
  ///@param Tw Temperature at the wall
  void mol_frac_to_part_dens(const RealVector& xb, CFreal press, CFreal Tw, RealVector& partrho);
  
  /// This function transform the molar fraction in partial density
  ///@param press  pressure at the wall
  ///@param yw     mass fractions at the wall
  ///@param dx     gradient of molar fractions normal to the wall
  ///@param Tw     temperature at the wall 
  ///@param flux   boundary flux (diffusion - mass production rate)
  void computeBoundaryFlux(CFreal& press,
			   RealVector& yw,
			   RealVector& dx,
			   CFreal& Tw,
			   RealVector& flux);
  
private: // data

  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter> m_inverter;
  
  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> m_numJacob;
  
  // the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;
  
  /// Gamma for N on the wall
  CFreal m_gammaN;
  
  /// Gamma for O on the wall
  CFreal m_gammaO;
  
  /// RealVector holding the input variables
  RealVector m_input;
  
  /// array to hold all vibrational temperatures and Te
  RealVector m_tVec;
  
  /// Massique fractions at the interior
  RealVector m_yi;
  
  /// Molar fractions at the interior
  RealVector m_xi;

  /// Mass fractions at the ghost
  RealVector m_yb;

  /// Mass fractions at the wall
  RealVector m_yw;
  
  /// perturbed Mass fractions at the wall
  RealVector m_yp;

  /// mass fractions at the ghost
  RealVector m_yg;

  /// Molar fractions at the ghost
  RealVector m_xb;
  
  /// Molar fractions at the wall
  RealVector m_xw;
 
  /// Vector of the molar fractions computed at the previous Newton step
  RealVector m_xwBkp;

  /// Difference between current and previous molar masses
  RealVector m_dxw;
 
  /// Molar fractions at the wall that we perturbe to compute 
  /// the diffusive flux
  RealVector m_xp;
  
  /// Vector of the molar masses
  RealVector m_mm;
  
  /// Vector of the derivative of the xb
  RealVector m_dx;
  
  /// Vector of the molar massestha we perturbe to compute 
  /// the diffusive flux
  RealVector m_dxp;

  /// Prodiction rate at the wall
  RealVector m_omega;
  RealVector m_omegap;

  /// Vector of the diffusive flux
  RealVector m_Diff_flux;
  
  /// Vector of the diffusive flux when we perturbe
  RealVector m_Diff_fluxp;
  
  /// Vector of the diffusive flux when we perturbe
  RealVector m_diff;
  
  /// Second hand right-hand side 
  RealVector m_b;
  
  /// Second hand right-hand side 
  RealVector m_zero;

  /// Vector of the impinging fluxes
  RealVector m_mcal;

  /// Partial densities
  RealVector m_partrho;
  
  /// temporary ghost rhoi
  RealVector m_rhoG;
  
  /// temporary gamma array
  RealVector m_gammaV;

  /// temporary array
  RealVector m_proj;

  /// Matrice of the diffusion coefficients (since we consider them vconstant here, it is not used)
  RealMatrix m_Dbin;
  
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
  
  // vector conatining the mu for all charged particles 
  std::vector<RealMatrix> m_muQ;
  
  /// Matrix of the derivative of the Diffusion flux
  /// minus the wall production wall
  RealMatrix m_a;
  
  /// Inverse matrix of the derivative of the Diffusion flux
  /// minus the wall production wall
  RealMatrix m_inva;
  
  /// Epsilon for numerical derivatives
  CFreal m_eps;

  /// Epsilon for controlling Newton iterations
  CFreal m_epsNewton; 
 
  /// Number of reactions catalysed by the wall
  CFuint m_nr;
  
  /// Number of Newton loop
  CFuint m_nl;
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;
  
  /// flag to choose Stefan-Maxwell diffusion fluxes
  bool m_useStefanMaxwell;
  
  /// a vector of reference values for molar fractions
  std::vector<CFreal> m_refXi;

}; // end of class NoSlipWallIsothermalNSrvtCatT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSrvtCatT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCatT_hh
