#ifndef COOLFluiD_Physics_NavierStokes_Euler2DLinearPuvt_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DLinearPuvt_hh

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
class Euler2DLinearPuvt : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DLinearPuvt(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DLinearPuvt();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler2DLinearPuvt

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DLinearPuvt_hh
