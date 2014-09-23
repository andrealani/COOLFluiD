#ifndef COOLFluiD_Numerics_FluctSplitNEQ_SplitStrategyNEQ_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_SplitStrategyNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace MathTools { class MatrixInverter; }
  
  namespace FluctSplitNEQ {

            
//////////////////////////////////////////////////////////////////////////////

/**
 * This class allows to bypass the function FluctuationSplitStrategy::computeConsistentStates()
 * for CNEQ flows
 *
 * @author Andrea Lani
 * @author Jesus Garicano Mena
 *
 */
template <class BASE>
class SplitStrategyNEQ : public BASE {
public:

  typedef Framework::MultiScalarTerm<COOLFluiD::Physics::NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SplitStrategyNEQ(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~SplitStrategyNEQ();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
protected: 
  
  /// Sets the current cell and calls the computation of the
  /// consistent state transformation.
  /// Notice that the parameters \chi_s and \kappa employed in this function, as defined in
  /// NASA CR 177489 (Liu & Vinokur), DIFFER from \beta, \phi and \gamma_s employed in
  /// NASA TP 2867(Gnoffo, Gupta and Shinn)
  virtual void setCurrentCell();

private: /// functions

void linearization_general_case();

void linearization_zero_gradP();

void eval_consistency( const CFreal& _beta_hat, const RealVector& _alpha_hat, const CFreal& _beta_bar, const RealVector& _alpha_bar );

private: /// data

  /// pointer to aerotermochemistry library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
  
  /// acquaintance of the update variables varset
  Common::SafePtr<Framework::ConvectiveVarSet> _updateVS;
  
  /// Vector to hold the mass fractions:
  RealVector _ys;

  /// Vector to hold species molecular masses:
  RealVector _mmasses;

  /// Vector to hold species gas constant (R/Mi):
  RealVector _RiGas;

  /// dpde at hat conditions
  CFreal _beta_hat;

  /// dpde at bar conditions
  CFreal _beta_bar;

  /// dpde at bar conditions
  CFreal _beta_bar_cons;

  RealVector _beta_hat_nodal;

  /// dpdRhoi at hat conditions
  RealVector _alpha_hat;

  /// dpdRhoi at bar conditions
  RealVector _alpha_bar;
  RealVector _alpha_bar_xx;
  RealVector _alpha_bar_yy;
  

  /// Translational-rotational energies of the species, at hat conditions
  RealVector _energyTr_hat;
  
  /// dpdRhoi at bar conditions
  RealVector _alpha_bar_cons;

  /// dpdRhoi at hat conditions
  RealVector _sqrtRho_atNodes;

  /// pressure gradient
  RealVector _gradP;

  /// gradients of sqrt(rho)*Y_s
  RealMatrix _gradSqrtRhoTimesYs;
  
  /// gradients of sqrt(rho)*Ux
  RealVector _gradSqrtRhoTimesU;
  
  /// gradients of sqrt(rho)*Vy
  RealVector _gradSqrtRhoTimesV;

  /// gradients of sqrt(rho)*H
  RealVector _gradSqrtRhoTimesH;
  
  /// gradients of sqrt(rho)*ev
  RealVector _gradSqrtRhoTimesEv;

//   /// gradient of rho
//   RealVector _gradRho;
// 
//   /// gradients of rho*Y_s
//   RealMatrix _gradRhoYs;
// 
//   /// gradients of rho*Ux
//   RealVector _gradRhoU;
// 
//   /// gradients of rho*Vy
//   RealVector _gradRhoV;
// 
//   /// gradients of rho*Ux*Ux
//   RealVector _gradRhoUSq;
// 
//   /// gradients of rho*Vy*Vy
//   RealVector _gradRhoVSq;
// 
//   /// gradients of rho*H
//   RealVector _gradRhoH;
// 
//   /// gradients of rho*ev
//   RealVector _gradRhoEv;

  /// Components of gradZ
  RealVector _dZdX;
  RealVector _dZdY;
  
  /// Components of gradJ
  RealVector _dJdX;
  RealVector _dJdY;

  /// gradients of rho*Y_s
  RealMatrix _consistent_gradRhoYs;

  /// gradients of rho*Y_s, projected
  RealVector _dRhoYsdL;

  /// gradients of rho*Ux*Ux
  RealVector _consistent_gradHalfRhoUSq;

  /// gradients of rho*H
  RealVector _consistent_gradRhoH;

  /// gradients of rho*ev
  RealVector _consistent_gradRhoEv;

  /// gradients of rho*E_tr
  RealVector _gradRhoE_tr;

  /// Vector to hold the averaged Z vector.
  RealVector _Zavg;

  /// Transformation matrix to compute dUdxi consistently
  RealMatrix _dPdZ;

  /// Transformation matrix to compute the gradients of rho_s, rhoU2, rhoH and rhoEv consistently
  RealMatrix _dJdZ;

  /// array storing the physical data
  RealVector _pData;

  /// array storing the physical data
  RealVector _pData0;
  /// array storing the physical data
  RealVector _pData1;
  /// array storing the physical data
  RealVector _pData2;

//   /// Needed for solving the 3x3 system of equations:
//   RealVector _rhs;
//   RealVector _solution;
//   RealMatrix _coefficients;
//   RealMatrix _inverseMat;
//   /// temporary data for holding the matrix inverter
//   MathTools::MatrixInverter* _inverter;

  /// Needed for solving the 2x2 system of equations:
  RealVector _rhs_Lagrange;
  RealVector _solution_Lagrange;
  RealMatrix _coefficients_Lagrange;
  RealMatrix _inverseMat_Lagrange;
  /// temporary data for holding the matrix inverter  
  MathTools::MatrixInverter* _inverter_Lagrange;

  /// for debugging purposes:
  RealVector _non_consistency;
  RealVector _consistency;

  std::string _message_state;
  std::string _message_consistency;
  std::string _message_betaBar_alphaBar;

  /// flag telling if debugging data will be stored
  CFuint _store_debugging_data;
    
}; // end of class SplitStrategyNEQ
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SplitStrategyNEQ.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_SplitStrategyNEQ_hh
