#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToCons_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToCons_hh

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
 *
 */
class Euler2DNEQRhoivtTvToCons : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQRhoivtTvToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~Euler2DNEQRhoivtTvToCons();
  
  /**
   * Transform a state into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
protected:
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
    
  /// Vector storing the elemental composition
  RealVector _ye;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;
 
  /// array with all different temperature
  RealVector _tDim;
  
  /// array with all different vibrational dimensional energies
  RealVector _evDim;

  /// array with partial densities 
  RealVector _rhoi; 

}; // end of class Euler2DNEQRhoivtTvToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToCons_hh
