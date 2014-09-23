#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingMFMaxwell2D_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingMFMaxwell2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeMultiFluidMHD/AUSMFluxMultiFluid.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM flux
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 *
 */
template <class UPDATEVAR>
class StegerWarmingMFMaxwell2D : public AUSMFluxMultiFluid<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  StegerWarmingMFMaxwell2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~StegerWarmingMFMaxwell2D();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data
   */
  virtual void setup();

protected:
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;

   /**
   * Compute the flux : implementation
   */
  virtual void computeMaxwellFlux();
  
   /**
   * Compute the Aplus Matrix
   */
  virtual void computeMatrixAplus();   
  
   /**
   * Compute the Aminus Matrix
   */
  virtual void computeMatrixAminus(); 
  
  
private:
  
  /// vector with the electromagnetic field variables (LEFT)
  RealVector _EMField_l;

  /// vector with the electromagnetic field variables (RIGHT)
  RealVector _EMField_r;
  
  /// A plus Matrix
  RealMatrix   _Aplus;
  
  /// A minus Matrix
  RealMatrix   _Aminus;
  
}; // end of class StegerWarmingMFMaxwell2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "StegerWarmingMFMaxwell2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingMFMaxwell2D_hh
