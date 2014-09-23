#ifndef COOLFluiD_Physics_MHD_MHD3DLinearPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DLinearPrim_hh

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
class MHD3DLinearPrim : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  MHD3DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~MHD3DLinearPrim();

   /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private:

  /// physical model
  Common::SafePtr<MHDTerm> _model;

}; // end of class MHD3DLinearPrim

//////////////////////////////////////////////////////////////////////////////

    } //  namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DLinearPrim_hh
