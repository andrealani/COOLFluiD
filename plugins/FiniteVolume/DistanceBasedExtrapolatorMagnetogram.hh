#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorMagnetogram_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorMagnetogram_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DistanceBasedExtrapolator.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object that uses a python
 * script to download and process a magnetogram
 *
 * @author Andrea Lani
 *
 */
class DistanceBasedExtrapolatorMagnetogram : public Framework::DistanceBasedExtrapolator<CellCenterFVMData> {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorMagnetogram(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorMagnetogram();

  /**
   * Set up private data needed by the computation
   */
  virtual void setup();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
protected:
  
  CFreal _sigma;

  CFreal _scaling_factor;

  /// path to the magnetogram data
  std::string _link;
  
  /// python command to use
  std::string _pyCommand;

  CFreal _Brefval;
  
  bool _runPyScript;
  
}; // end of class DistanceBasedExtrapolatorMagnetogram

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorMagnetogram_hh
