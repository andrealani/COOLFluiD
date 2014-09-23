#ifndef COOLFluiD_Physics_NavierStokes_Euler3DLinearCons_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DLinearCons_hh

//////////////////////////////////////////////////////////////////////////////

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
class Euler3DLinearCons : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler3DLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler3DLinearCons();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler3DLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NavierStokes


  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DLinearCons_hh
