#ifndef COOLFluiD_Physics_LTE_Euler2DLinearPrvtLTE_hh
#define COOLFluiD_Physics_LTE_Euler2DLinearPrvtLTE_hh

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
     }
    
    namespace LTE {
      
//////////////////////////////////////////////////////////////////////////////

 /**
  * This class computes linearized values accepting [p rhoU rhoV T]
  * variables and setting the EulerData corresponding to a Euler
  * with LTE model.
  *
  * @author Andrea Lani
  */
class Euler2DLinearPrvtLTE : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DLinearPrvtLTE(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DLinearPrvtLTE();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private:

 /// acquaintance of the PhysicalModel
  Common::SafePtr<NavierStokes::EulerTerm> _model;
  
  /// array with density, enthalpy and energy
  RealVector _dhe;

}; // end of class Euler2DLinearPrvtLTE

//////////////////////////////////////////////////////////////////////////////

    } //  namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_Euler2DLinearPrvtLTE_hh
