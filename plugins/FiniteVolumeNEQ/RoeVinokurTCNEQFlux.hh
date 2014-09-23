#ifndef COOLFluiD_Numerics_FiniteVolume_RoeVinokurTCNEQFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeVinokurTCNEQFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux for thermo-chemical non-equilibrium applications
 *
 * @author Andrea Lani
 *
 */
template <class BASE , class UPDATEVAR>
class RoeVinokurTCNEQFlux : public BASE {
public:
  
  typedef Framework::MultiScalarTerm<COOLFluiD::Physics::NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Constructor
   */
  RoeVinokurTCNEQFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeVinokurTCNEQFlux();

  /**
   * Set up private data
   */
  virtual void setup();
    
private:
  
  /**
   * Perform ad-hoc linearization
   */
  void linearize();
  
private:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> _upVar;
  
  /// physical model data
  RealVector _lData;

  /// physical model data
  RealVector _rData;
  
//   /// species molar masses
//   RealVector _mmasses;
   
  /// f_i coefficients
  RealVector _fcoeff;
  
  /// total variation of partial densities over corresponding molar masses 
  RealVector _dRhoiOvMM;

//   /// pointer to aerotermochemistry library
//   Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
// 
//   /// acquaintance of the model
//   Common::SafePtr<NEQTerm> _model;
// 
//   /// acquaintance of the update variables varset
//   Common::SafePtr<Framework::ConvectiveVarSet> _updateVS;
// 
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

  RealVector _beta_hat_nodal;

  /// dpdRhoi at hat conditions
  RealVector _alpha_hat;

  /// dpdRhoi at bar conditions
  RealVector _alpha_bar;

  /// Translational-rotational energies of the species, at Left and Right states
  RealVector _energyTr_L;
  RealVector _energyTr_R;

  /// Translational-rotational energies of the species, at hat conditions
  RealVector _energyTr_hat;
  
  /// dpdRhoi at bar conditions
  RealVector _alpha_bar_cons;

  /// dpdRhoi at hat conditions
  RealVector _sqrtRho_atNodes;

  /// pressure gradient
  CFreal _DeltaP;

  /// gradients of sqrt(rho)*Y_s
  RealVector _DeltaSqrtRhoTimesYs;

  /// gradients of sqrt(rho)*Ux
  CFreal _DeltaSqrtRhoTimesU;

  /// gradients of sqrt(rho)*Vy
  CFreal _DeltaSqrtRhoTimesV;

  /// gradients of sqrt(rho)*H
  CFreal _DeltaSqrtRhoTimesH;

  /// gradients of sqrt(rho)*ev
  CFreal _DeltaSqrtRhoTimesEv;

  /// Components of gradZ
  RealVector _DeltaZ;
  
  /// Components of gradJ
  RealVector _DeltaJ;  

  /// gradients of rho*Y_s
  RealVector _consistent_DeltaRhoYs;

  /// gradient of rho*(Ux*Ux + Uy*Uy)
  CFreal _consistent_DeltaHalfRhoUSq;

  /// gradients of rho*H
  CFreal _consistent_DeltaRhoH;

  /// gradients of rho*ev
  CFreal _consistent_DeltaRhoEv;

  /// Vector to hold the averaged Z vector.
  RealVector _ZRight;
  RealVector _ZLeft;

  /// Vector to hold the averaged Z vector.
  RealVector _Zavg;
  
  /// Transformation matrix to compute the gradients of rho_s, rhoU2, rhoH and rhoEv consistently
  RealMatrix _dJdZ;

  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
  /// pointer to aerotermochemistry library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

        
}; // end of class RoeVinokurTCNEQFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "RoeVinokurTCNEQFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeVinokurTCNEQFlux_hh
