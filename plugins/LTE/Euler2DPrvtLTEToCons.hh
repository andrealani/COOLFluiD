#ifndef COOLFluiD_Physics_LTE_Euler2DPrvtLTEToCons_hh
#define COOLFluiD_Physics_LTE_Euler2DPrvtLTEToCons_hh

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
class Euler2DPrvtLTEToCons : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DPrvtLTEToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DPrvtLTEToCons();

  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<NavierStokes::EulerTerm> _model;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
}; // end of class Euler2DPrvtLTEToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_Euler2DPrvtLTEToCons_hh
