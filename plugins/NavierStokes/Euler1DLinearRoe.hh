#ifndef COOLFluiD_Physics_NavierStokes_Euler1DLinearRoe_hh
#define COOLFluiD_Physics_NavierStokes_Euler1DLinearRoe_hh

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
  * @author Alessandro Munafo'
  */
class Euler1DLinearRoe : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler1DLinearRoe(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler1DLinearRoe();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: // data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler1DLinearRoe

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler1DLinearRoe_hh
