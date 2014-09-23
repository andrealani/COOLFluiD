#ifndef COOLFluiD_Physics_NavierStokes_Euler2DLinearRhovt_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DLinearRhovt_hh

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
class Euler2DLinearRhovt : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DLinearRhovt(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DLinearRhovt();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler2DLinearRhovt

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DLinearRhovt_hh
