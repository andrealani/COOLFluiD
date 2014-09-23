#ifndef COOLFluiD_Physics_MHD_MHD2DProjectionLinearCons_hh
#define COOLFluiD_Physics_MHD_MHD2DProjectionLinearCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace MHD {
      class MHDProjectionTerm;

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Andrea Lani
  * @author Mehmet Sarp Yalim
  */
class MHD2DProjectionLinearCons : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  MHD2DProjectionLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~MHD2DProjectionLinearCons();

   /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private:

  /// physical model
  Common::SafePtr<MHDProjectionTerm> _model;

}; // end of class MHD2DProjectionLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DProjectionLinearCons_hh
