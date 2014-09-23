#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler1DTtPtYiMultiTvTe_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler1DTtPtYiMultiTvTe_hh

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
  * This class implements a subsonic inlet in terms of tTotal, pTotal, Yi, the various Tv and Te
  *
  * @author Alessandro Munafò
  *
  */
   
class SubInletEuler1DTtPtYiMultiTvTe : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler1DTtPtYiMultiTvTe(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletEuler1DTtPtYiMultiTvTe();

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

  /// molar masses
  RealVector _mm;

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;
 
  /// Yi mass fraction
  std::vector<CFreal> _Yi;

  /// Vibrational temperatures
  std::vector<CFreal>  _Tv;

  /// Electron temperature
  CFreal     _Te;
  
}; // end of class SubInletEuler1DTtPtYiTv

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler1DTtPtYiTv_hh
