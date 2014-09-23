#ifndef COOLFluiD_Physics_LTE_IncompEuler2DdPuvtLTEToCons_hh
#define COOLFluiD_Physics_LTE_IncompEuler2DdPuvtLTEToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace NavierStokes {
      class IncompEulerTerm;
    }
    
    namespace LTE {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [p2 u v T] to conservative variables
 *
 * @author Radek Honzattko
 *
 */
class IncompEuler2DdPuvtLTEToCons : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  IncompEuler2DdPuvtLTEToCons
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~IncompEuler2DdPuvtLTEToCons();
  
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

  /// acquaintance of the PhysicalModel
  Common::SafePtr<NavierStokes::IncompEulerTerm> _model;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
}; // end of class IncompEuler2DdPuvtLTEToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_IncompEuler2DdPuvtLTEToCons_hh
