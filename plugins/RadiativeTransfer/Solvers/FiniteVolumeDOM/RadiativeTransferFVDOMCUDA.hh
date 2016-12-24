#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferFVDOMCUDA_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferFVDOMCUDA_hh

///////////////////////////////////////////////////////////////////////////

#include "RadiativeTransfer/Solvers/FiniteVolumeDOM/RadiativeTransferFVDOM.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace RadiativeTransfer {
      
//////////////////////////////////////////////////////////////////////////

/**
 * This class compute the radiative heat transfer using a Finite Volume algorithm
 *
 * @author Andrea Lani
 */
class RadiativeTransferFVDOMCUDA : public RadiativeTransferFVDOM {
public:

  /**
   * Constructor
   */
  RadiativeTransferFVDOMCUDA(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RadiativeTransferFVDOMCUDA();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);  

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();
  
private: //function
  
  /// Compute radiative fluxes by looping over bins
  virtual void loopOverBins(const CFuint startBin, 
			    const CFuint endBin, 
			    const CFuint startDir,
			    const CFuint endDir);
  
  /// Compute radiative fluxes by looping over directions
  virtual void loopOverDirs(const CFuint startBin, 
			    const CFuint endBin, 
			    const CFuint startDir,
			    const CFuint endDir);
  
private: //data
  
  /// face-cell connectivity
  Framework::LocalArray<CFint>::TYPE m_faceCell;

  /// number of faces per cell
  Framework::LocalArray<CFuint>::TYPE m_nbFacesInCell;
  
  /// Exponent for the radiation of oppacity table
  Framework::LocalArray<CFreal>::TYPE m_InDir;
  
  /// qx for each direction
  Framework::LocalArray<CFreal>::TYPE m_qxDir;
  
  /// qz for each direction
  Framework::LocalArray<CFreal>::TYPE m_qyDir;
  
  /// qz for each direction
  Framework::LocalArray<CFreal>::TYPE m_qzDir;
  
  /// divq for each direction
  Framework::LocalArray<CFreal>::TYPE m_divqDir;
  
  /// name of the algorithm to use for computing Q
  std::string m_qAlgoName;
  
}; // end of class RadiativeTransferFVDOMCUDA
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_RadiativeTransferFVDOMCUDA_hh



