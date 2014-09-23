#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQLinearRhoivtCNEQ_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQLinearRhoivtCNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
    template <class TERM> class MultiScalarTerm;
  }

  namespace Physics {

    namespace NavierStokes {
      class EulerTerm;
    }
    
    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs a linearization of phyical quantities in the case of 
  * thermo-chemical NEQ
  *
  * @author Andrea Lani
  * @author Jesus Garicano Mena
  * 
  */
class Euler2DNEQLinearRhoivtCNEQ : public Framework::JacobianLinearizer {
public:
    
  /**
   * Default constructor without arguments
   */
  Euler2DNEQLinearRhoivtCNEQ(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQLinearRhoivtCNEQ();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: //data
  
  /// acquaintance of the PhysicalModel
  Common::SafePtr<Framework::MultiScalarTerm<NavierStokes::EulerTerm> > _model;
    
  /// Vector storing the elemental composition
  RealVector _ye;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
}; // end of class Euler2DNEQLinearRhoivtCNEQ

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQLinearRhoivtCNEQ_hh
