#ifndef COOLFluiD_Physics_GReKO_Euler2DGReKOLinearRoe_hh
#define COOLFluiD_Physics_GReKO_Euler2DGReKOLinearRoe_hh

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

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Khalil Bensassi 
  */

class Euler2DGReKOLinearRoe : public Framework::JacobianLinearizer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerGReKOTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DGReKOLinearRoe(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DGReKOLinearRoe();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: // data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerGReKOTerm> _model;

}; // end of class Euler2DGReKOLinearRoe

//////////////////////////////////////////////////////////////////////////////

    } //  namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_Euler2DGReKOLinearRoe_hh
