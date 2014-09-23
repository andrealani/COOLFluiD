#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQLinearCons_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQLinearCons_hh

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
  */
class Euler2DNEQLinearCons : public Framework::JacobianLinearizer {
public:
    
  /**
   * Default constructor without arguments
   */
  Euler2DNEQLinearCons(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQLinearCons();

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
  
  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;
  
  /// array with all different vibrational dimensional energies
  RealVector _evDim;
  
 /// species molar masses
  RealVector _mmasses;
    
}; // end of class Euler2DNEQLinearCons

//////////////////////////////////////////////////////////////////////////////

    } //  namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQLinearCons_hh
