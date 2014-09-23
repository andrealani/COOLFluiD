#ifndef COOLFluiD_Numerics_LESDataProcessing_GradientComputerFVMCC_hh
#define COOLFluiD_Numerics_LESDataProcessing_GradientComputerFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "GradientComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

/** 
 * This class calculates gradients in primitive variables using the 
 * Green-Gauss theorem for cell-centres.
 * 
 * @author Willem Deconinck
 */
class GradientComputerFVMCC : public GradientComputer {

  
public: // functions

  /** 
   * Defines the Config Options of this class
   * @param    options   an OptionList where to add the options
   */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  GradientComputerFVMCC(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~GradientComputerFVMCC();

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /** 
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "GradientComputerFVMCC";
  }
  
  virtual void compute(std::vector<RealVector>& gradients, const CFuint& iCell) ;
  
  virtual CFreal getVolume(const CFuint& iCell);
  
  virtual CFreal getVolumeAdim(const CFuint& iCell);
  
protected: // helper functions

private: // data

  // Connectivity of cells to cell faces
  Common::SafePtr<Common::ConnectivityTable<CFuint> > m_cellFaces;

  // To find which TRS the face belongs to
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs;

  // To build faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> > m_geoBuilder;

  RealVector m_primState;

  /// current face
  Framework::GeometricEntity* m_currFace;

  // reference volume and reference area
  CFreal m_refVol;
  CFreal m_refArea;

  // temporary averageState storage
  RealVector m_avState;

  // temporary normal storage
  RealVector m_normal;
  
protected: // data



}; // end of class GradientComputerFVMCC

//////////////////////////////////////////////////////////////////////////////

  	} // namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_GradientComputerFVMCC_hh
