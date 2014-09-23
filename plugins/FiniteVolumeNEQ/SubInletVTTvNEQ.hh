#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletVTTvNEQ_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletVTTvNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

 /**
   * This class represents a command that applies the Subsonic Inlet BC
   * in case of update variables [rho_i v T Tv]
   *  
   * @author Andrea Lani
   *
   */
template <class MODEL>
class SubInletVTTvNEQ : public FiniteVolume::FVMCC_BC {

public: 
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SubInletVTTvNEQ(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~SubInletVTTvNEQ();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
private:
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// pointer to the convective term
  Common::SafePtr<MODEL> m_model;
  
  /// mass species fractions in the internal cell
  RealVector m_ysIn;
  
  /// physical data array
  RealVector m_pData;
  
  /// inlet values to be fixed [V T Tv]
  /// the partial pressures will be extrapolated from inside 
  std::vector<CFreal> m_vTTv;
  
}; // end of class NoSlipWallIsothermalNSPvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubInletVTTvNEQ.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletVTTvNEQ_hh
