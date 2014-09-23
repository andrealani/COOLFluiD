#ifndef COOLFluiD_Physics_NavierStokes_Euler2DLinearRoe_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DLinearRoe_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/BadValueException.hh"
#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace NavierStokes {
      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Andrea Lani
  */
class Euler2DLinearRoe : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DLinearRoe(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DLinearRoe();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: // data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler2DLinearRoe

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DLinearRoe_hh
