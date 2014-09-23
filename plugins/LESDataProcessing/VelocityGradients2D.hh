#ifndef COOLFluiD_Numerics_LESDataProcessing_VelocityGradients2D_hh
#define COOLFluiD_Numerics_LESDataProcessing_VelocityGradients2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "TurbulenceFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {
		  
		  class LESProcessingData;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This strategy class computes the Velocity Gradients in 2D
 * and stores it in the socket named "VelocityGradients"
 * @author Willem Deconinck
 */
class VelocityGradients2D : public TurbulenceFunction {

public: // functions
  
  enum {RHO,U,V,P};
  
	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  VelocityGradients2D(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~VelocityGradients2D();

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
    return "VelocityGradients2D";
  }
  
  virtual void setSocketName()
  {
    m_socketName = "VelocityGradients";
  }

  virtual void compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell);

protected:
  
  /// Which velocity components to store in the socket
  std::vector<std::string> m_velocityComponents;
  
  /// The velocity components in index form
  std::vector<CFuint> m_comp;
  
  /// The size of m_velocityComponents and m_comp
  CFuint m_nbVelocityComponents;
  
  /// The number of scalars to be stored per cell in the socket
  CFuint m_stride;
  
  /// The dimension.
  CFuint m_dim;

}; // end of class VelocityGradients2D

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_VelocityGradients2D_hh
