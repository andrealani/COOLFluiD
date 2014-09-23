	#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DTtPtAlphaYiMultiTv_hh
	#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DTtPtAlphaYiMultiTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  

  namespace Physics {
    namespace NavierStokes {
      template <class BASEVS> class MultiScalarVarSet;
      class Euler2DVarSet;
    }
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
 /**
  * This class implements a subsonic inlet given tTotal, pTotal, alpha, Yi, Tv
  *
  * @author Alessandro Munafo'
  * @author Andrea Lani
  *
  *
  */
   
class SubInletEuler2DTtPtAlphaYiMultiTv : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler2DTtPtAlphaYiMultiTv(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletEuler2DTtPtAlphaYiMultiTv();

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
  Common::SafePtr<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DVarSet> > _varSet;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// physical model data
  RealVector _dataInnerState;

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;

  /// alpha
  CFreal     _alpha;
  
  /// vibrational temperatures
  std::vector<CFreal>  _Tv;
  
  /// Yi mass fraction
  std::vector<CFreal> _Yi;
 
}; // end of class SubInletEuler2DTtPtAlphaYiMultiTv

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DTtPtAlphaYiMultiTv_hh
