#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the Isothermal wall
   * BC in case of update variables including velocity components and temperature
   * [p v T y_i] or [rho_i v T] etc.
   *  
   * @author Andrea Lani
   *
   */
template <class MODEL>
class NoSlipWallIsothermalNSrvt : public NoSlipWallIsothermalNSvt<MODEL> {

public: 

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSrvt(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NoSlipWallIsothermalNSrvt();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
protected:
  
  /**
   * Apply boundary condition on the given ghost state
   * @param innerState state vector associated to the internal cell
   * @param ghostState ghost state vector associated to the ghost cell
   * @param 
   */
  virtual void setGhostStateImpl(const Framework::State& innerState, 
				 Framework::State& ghostState);
  
protected:
  
  /// physico-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// flag telling if the state has partial densities
  bool m_stateHasPartialDensities;
  
  /// number of species
  CFuint m_nbSpecies;
  
  /// number of vibrational temperatures
  CFuint m_nbTv;
  
  /// roto-translational and vibrational temperatures in the ghost state
  RealVector m_ghostTTvib;
  
  /// roto-translational and vibrational temperatures in the inner state
  RealVector m_innerTTvib;
   
}; // end of class NoSlipWallIsothermalNSPvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSrvt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvt_hh
