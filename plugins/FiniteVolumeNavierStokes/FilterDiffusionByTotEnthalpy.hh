#ifndef COOLFluiD_Numerics_FiniteVolume_FilterDiffusionByTotEnthalpy_hh
#define COOLFluiD_Numerics_FiniteVolume_FilterDiffusionByTotEnthalpy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_EquationFilter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a generic FVMCC equation filter for diffusion terms
 *
 * @author Andrea Lani
 *
 */
class FilterDiffusionByTotEnthalpy : public FVMCC_EquationFilter {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FilterDiffusionByTotEnthalpy(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FilterDiffusionByTotEnthalpy();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Unsetup up private data to prepare the simulation
   */
  virtual void unsetup();
   
  /// Reset all data before starting a new iteration
  virtual void reset();
  
  /// Filter the equation on the current geometric entity
  virtual bool filterOnGeo(Framework::GeometricEntity *const geo);
  
protected:
  
  /// factor
  CFreal m_factorH;
  
  /// free stream total enthalpy
  CFreal m_freeStreamH;
  
  /// deviation from free stream total enthalpy
  CFreal m_deviationH;
  
  /// flag to telll if to interpret data coming from postprocessing
  bool m_fromPostProcessing;
  
}; // end of class FilterDiffusionByTotEnthalpy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FilterDiffusionByTotEnthalpy_hh
