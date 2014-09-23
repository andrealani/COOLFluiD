#ifndef COOLFluiD_Numerics_LESDataProcessing_Vorticity3D_hh
#define COOLFluiD_Numerics_LESDataProcessing_Vorticity3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "TurbulenceFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {
		  
		  class LESProcessingData;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This strategy class computes the vorticity in 3D
 * and stores it in the socket named "Vorticity"
 * 
 * @author Willem Deconinck
 */
class Vorticity3D : public TurbulenceFunction {

public: // functions
  
  enum {RHO,U,V,W,P};
  
	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  Vorticity3D(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~Vorticity3D();

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
    return "Vorticity3D";
  }
  
  virtual void setSocketName()
  {
    m_socketName = "Vorticity";
  }

  virtual void compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell);

protected:

}; // end of class Vorticity3D

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_Vorticity3D_hh
