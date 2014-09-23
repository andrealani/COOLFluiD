#ifndef COOLFluiD_Physics_NEQ_Euler1DNEQRhoivtTvToCons_hh
#define COOLFluiD_Physics_NEQ_Euler1DNEQRhoivtTvToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [rhoi u v T Tvi] to conservative variables
 *
 * @author Andrea Lani
 * @author Alessandro Munafò
 *
 */
class Euler1DNEQRhoivtTvToCons : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler1DNEQRhoivtTvToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler1DNEQRhoivtTvToCons();
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
    
  /// Vector storing the elemental composition
  RealVector _ye;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;
  
  /// array with all different vibrational dimensional energies
  RealVector _evDim;
  
}; // end of class Euler1DNEQRhoivtTvToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler1DNEQRhoivtTvToCons_hh
