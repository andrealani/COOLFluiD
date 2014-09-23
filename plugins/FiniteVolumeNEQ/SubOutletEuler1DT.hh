#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEuler1DT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEuler1DT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  

  namespace Physics {
    namespace NavierStokes {
      template <class BASEVS> class MultiScalarVarSet;
      class Euler1DVarSet;
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
  * imposing the temperature
  *
  * @author Alessandro Munafò
  *
  */
   
class SubOutletEuler1DT : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletEuler1DT(const std::string& name);

  /**
   * Default destructor
   */
  ~SubOutletEuler1DT();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private: // data
  
  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler1DVarSet> > _varSet;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// physical model data
  RealVector _dataInnerState;
  
  /// total pressure
  CFreal     _tOut;
  
}; // end of class SubOutletEuler1DT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEuler1DT_hh
