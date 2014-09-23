#ifndef COOLFluiD_Physics_LTE_Euler2DRhovtXiToPuvtLTE_hh
#define COOLFluiD_Physics_LTE_Euler2DRhovtXiToPuvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace NavierStokes {
      class EulerTerm;
    }
    
    namespace LTE {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from [rho v T Xi] to 
 * [p u v T] variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DRhovtXiToPuvtLTE : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DRhovtXiToPuvtLTE(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DRhovtXiToPuvtLTE();

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

  /// acquaintance of the PhysicalModel
  Common::SafePtr<NavierStokes::EulerTerm> _model;
  
  /// array storing density enthalpy energy
  RealVector _dhe;

  /// Vector storing the mass fractions
  RealVector _ye;
  
  /// Vector storing the molar fractions
  RealVector _xe;

  /// Vector storing the molar masses
  RealVector _mmasses;

}; // end of class Euler2DRhovtXiToPuvtLTE

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_Euler2DRhovtXiToPuvtLTE_hh
