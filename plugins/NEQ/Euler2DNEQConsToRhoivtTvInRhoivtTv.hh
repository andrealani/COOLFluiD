#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQConsToRhoivtTvInRhoivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQConsToRhoivtTvInRhoivtTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * Transformer from conservative to [rho_i, u, v, T, Tv] variables, with the
 * matrix evaluated in [rho_i, u, v, T, Tv].
 *
 * This is the inverse of Euler2DNEQRhoivtTvToConsInRhoivtTv and is what
 * SolToUpdateInUpdateMatTrans needs when the solution variables are
 * conservative and the update variables are RhoivtTv.
 *
 * The matrix is built analytically by inverting the forward one rather than
 * by numerical inversion, so it stays exact and costs nothing extra.
 *
 * NOT usable with PLATO as the property library. Like the forward transformer
 * it needs dEdT, dEvTv, dRhoEdRhoi and dRhoEvdRhoi from
 * PhysicalChemicalLibrary::ExtraData, and PlatoLibrary::setDensityEnthalpyEnergy
 * throws NotImplementedException when asked to store them. It is also not on
 * the critical path for the FR TCNEQ setup: the only consumer of
 * SolToUpdateInUpdateMatTrans is FinalizeRHS, and every testcase in this repo
 * runs with FinalizeRHSCom = Null.
 *
 * @author Rayan Dhib
 */
class Euler2DNEQConsToRhoivtTvInRhoivtTv : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DNEQConsToRhoivtTvInRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQConsToRhoivtTvInRhoivtTv();

  /**
   * Set the transformation matrix from a given state
   */
  void setMatrix(const RealVector& state);

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data

  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;

  /// Vector storing the elemental composition
  RealVector _ys;

  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;

  /// array with all different vibrational dimensional energies
  RealVector _evDim;

  /// array to store density, enthalpy and energy
  RealVector _dhe;

}; // end of class Euler2DNEQConsToRhoivtTvInRhoivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQConsToRhoivtTvInRhoivtTv_hh
