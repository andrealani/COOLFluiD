#ifndef COOLFluiD_Physics_NavierStokes_EulerPvtToCons_hh
#define COOLFluiD_Physics_NavierStokes_EulerPvtToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [p v T] to conservative variables
 *
 * @author Andrea Lani
 *
 */
class EulerPvtToCons : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  EulerPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~EulerPvtToCons();
  
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
  Common::SafePtr<EulerTerm> _model;
  
  /// density
  CFreal m_rho;
  
}; // end of class EulerPvtToCons
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_EulerPvtToCons_hh
