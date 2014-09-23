#ifndef COOLFluiD_Physics_KOmega_Euler2DKOmegaLinearCons_hh
#define COOLFluiD_Physics_KOmega_Euler2DKOmegaLinearCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/BadValueException.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Lilla Koloszar
  */

class Euler2DKOmegaLinearCons : public Framework::JacobianLinearizer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerKOmegaTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DKOmegaLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DKOmegaLinearCons();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: // data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerKOmegaTerm> _model;

}; // end of class Euler2DKOmegaLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_Euler2DKOmegaLinearCons_hh
