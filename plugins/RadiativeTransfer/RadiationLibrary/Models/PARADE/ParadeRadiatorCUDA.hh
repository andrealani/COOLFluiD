#ifndef COOLFluiD_RadiativeTransfer_ParadeRadiatorCUDA_hh
#define COOLFluiD_RadiativeTransfer_ParadeRadiatorCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "RadiativeTransfer/RadiationLibrary/Models/PARADE/ParadeRadiator.hh"

///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a parallel radiator using the PARADE library and 
 * taking advantage of heterogeneous computing (CPU/GPU) for implementing 
 * the most time consuming functions 
 *
 * @author Andrea Lani
 *
 */
class ParadeRadiatorCUDA : public ParadeRadiator {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  ParadeRadiatorCUDA(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ParadeRadiatorCUDA();
    
  /**
   * Setups the data of the library
   */
  virtual void setup();
  
  /**
   * Unsetups the data of the library
   */
  virtual void unsetup();
  
protected:
  
  /// apply the binning method to reduce the spectral data
  virtual void computeBinning();
  
  /// apply the banding method to reduce the spectral data
  virtual void computeBanding();
  
  /// apply the binning/banding method to reduce the spectral data
  virtual void computeBinningBanding();
  
  /// compute the averaged bins/bands
  virtual void computeAveragedBins(const CFuint nbBinsre, 
				   const CFuint testID,
				   Framework::LocalArray<CFreal>::TYPE& vctBins);
  
}; // end of class ParadeRadiatorCUDA
    
//////////////////////////////////////////////////////////////////////////////
    
  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_ParadeRadiatorCUDA_hh
