#ifndef COOLFluiD_Physics_SA_Euler2DSALinearCons_hh
#define COOLFluiD_Physics_SA_Euler2DSALinearCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/BadValueException.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {
    
    namespace SA {
      
//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Andrea Lani
  */
class Euler2DSALinearCons : public Framework::JacobianLinearizer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerSATerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DSALinearCons(Common::SafePtr<Framework::PhysicalModel> model);
  
  /**
   * Default destructor
   */
  ~Euler2DSALinearCons();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: // data
  
  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerSATerm> _model;
  
}; // end of class Euler2DSALinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_Euler2DSALinearCons_hh
