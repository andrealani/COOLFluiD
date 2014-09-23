#ifndef COOLFluiD_Physics_LTE_Euler3DPvtLTEToCons_hh
#define COOLFluiD_Physics_LTE_Euler3DPvtLTEToCons_hh

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
 * This class represents a transformer of variables from primitive
 * [p u v T] to conservative variables
 *
 * @author Andrea Lani
 *
 */
class Euler3DPvtLTEToCons : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler3DPvtLTEToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~Euler3DPvtLTEToCons();
  
  /**
   * Transform a state into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<NavierStokes::EulerTerm> _model;
  
  /// array storing density enthalpy energy
  RealVector _dhe;
  
}; // end of class Euler3DPvtLTEToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_Euler3DPvtLTEToCons_hh
