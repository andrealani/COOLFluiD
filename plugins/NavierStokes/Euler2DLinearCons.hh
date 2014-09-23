#ifndef COOLFluiD_Physics_NavierStokes_Euler2DLinearCons_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DLinearCons_hh

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
class Euler2DLinearCons : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DLinearCons();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler2DLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DLinearCons_hh
