#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEulerP_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEulerP_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  

  namespace Physics {
    namespace NavierStokes {
      template <class BASEVS> class MultiScalarVarSet;
    }
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
 /**
  * This class implements a subsonic outlet boundary conditions by 
  * imposing the pressure
  *
  * @author Alessandro Munafò
  * @author Andrea Lani
  *
  */
template <typename UPDATEVAR>      
class SubOutletEulerP : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletEulerP(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubOutletEulerP();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
protected: // data
  
  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::MultiScalarVarSet<UPDATEVAR> > _varSet;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// flag telling if the state has partial densities
  bool _stateHasPartialDensities;
  
  /// physical model data
  RealVector _dataInnerState;
  
  /// molar masses
  RealVector _mm;
    
  /// R*T*y_i/M_i per each species
  RealVector _RTYovM;
  
  /// total pressure
  CFreal     _pOut;
  
}; // end of class SubOutletEulerP

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubOutletEulerP.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEulerP_hh
