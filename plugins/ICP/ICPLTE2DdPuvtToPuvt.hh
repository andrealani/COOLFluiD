#ifndef COOLFluiD_Physics_ICP_ICPLTE2DdPuvtToPuvt_hh
#define COOLFluiD_Physics_ICP_ICPLTE2DdPuvtToPuvt_hh

//////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ICP {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [p2 u v T] to conservative variables
 *
 * @author Radek Honzattko
 *
 */
class ICPLTE2DdPuvtToPuvt : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
  
  /**
   * Default constructor without arguments
   */
  ICPLTE2DdPuvtToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~ICPLTE2DdPuvtToPuvt();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<PTERM> _icpModel;
  
  /// array storing density enthalpy energy
  RealVector _dhe;
  
}; // end of class ICPLTE2DdPuvtToPuvt

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPLTE2DdPuvtToPuvt_hh
