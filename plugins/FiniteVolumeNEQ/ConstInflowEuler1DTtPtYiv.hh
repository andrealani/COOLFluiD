#ifndef COOLFluiD_Numerics_FiniteVolume_ConstInflowEuler1DTtPtYiv_hh
#define COOLFluiD_Numerics_FiniteVolume_ConstInflowEuler1DTtPtYiv_hh

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
  * This class implements a constant mass flow given tTotal, pTotal, Yi and u
  *
  * @author Alessandro Munafò
  *
  */
   
class ConstInflowEuler1DTtPtYiv : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  ConstInflowEuler1DTtPtYiv(const std::string& name);

  /**
   * Default destructor
   */
  ~ConstInflowEuler1DTtPtYiv();

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

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;

  // velocity 
  CFreal    _u;
  
  /// Yi mass fraction
  std::vector<CFreal> _Yi;
  
}; // end of class SubInletEuler1DTtPtYi

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ConstInflowEuler1DTtPtYiv_hh
