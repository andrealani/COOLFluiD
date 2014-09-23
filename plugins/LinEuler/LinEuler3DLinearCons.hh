#ifndef COOLFluiD_Physics_LinEuler_LinEuler3DLinearCons_hh
#define COOLFluiD_Physics_LinEuler_LinEuler3DLinearCons_hh

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
class LinEuler3DLinearCons : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  LinEuler3DLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~LinEuler3DLinearCons();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);


private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<LinEulerTerm> _model;

}; // end of class LinEuler3DLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Euler_Euler2DLinearCons_hh
