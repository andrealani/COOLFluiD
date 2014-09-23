#ifndef COOLFluiD_Physics_MHD_MHD2DLinearCons_hh
#define COOLFluiD_Physics_MHD_MHD2DLinearCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace MHD {
      class MHDTerm;

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Andrea Lani
  */
class MHD2DLinearCons : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  MHD2DLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~MHD2DLinearCons();

   /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private:

  /// physical model
  Common::SafePtr<MHDTerm> _model;

}; // end of class MHD2DLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DLinearCons_hh
