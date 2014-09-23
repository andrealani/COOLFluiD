/* CHANGED FOR UNHOMOGENEOUS MEAN FLOW */

#ifndef COOLFluiD_Physics_LinEuler_LinEuler2DLinearCons_hh
#define COOLFluiD_Physics_LinEuler_LinEuler2DLinearCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;
  }

  namespace Physics {

    namespace LinearizedEuler {

      class LinEulerTerm;

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Lilla Koloszar
  * @author Nadege Villedieu
  * @author Tomas Kopacek
  */
class LinEuler2DLinearCons : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  LinEuler2DLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~LinEuler2DLinearCons();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);


private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<LinEulerTerm> _model;

}; // end of class LinEuler2DLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Euler_Euler2DLinearCons_hh
