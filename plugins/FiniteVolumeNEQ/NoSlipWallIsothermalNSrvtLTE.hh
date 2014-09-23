#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtLTE_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"

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
   * [p v T y_i] or [rho_i v T] etc. Local Thermodynamic equilibrium is assumed 
   * at the wall, meaning that the species composition and the density at the wall
   * must be the equilibrium ones
   *  
   * @author Andrea Lani
   *
   */
template <class MODEL>
class NoSlipWallIsothermalNSrvtLTE : public NoSlipWallIsothermalNSrvt<MODEL> {

public: 
 
  /**
   * Constructor
   */
  NoSlipWallIsothermalNSrvtLTE(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NoSlipWallIsothermalNSrvtLTE();
  
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

private:
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// pointer to the convective term
  Common::SafePtr<MODEL> m_model;
  
  /// electron-vibrational temperatures
  RealVector m_tVec;
  
  /// mass species fractions in the internal cell
  RealVector m_ysIn;
  
  /// mass species fractions at equilibrium for local (p,T,Xe)
  RealVector m_ysEq;
  
  /// array of species molar fractions
  RealVector m_xs;

}; // end of class NoSlipWallIsothermalNSPvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSrvtLTE.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtLTE_hh
