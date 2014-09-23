#ifndef COOLFluiD_Numerics_LESDataProcessing_Qcriterion3D_hh
#define COOLFluiD_Numerics_LESDataProcessing_Qcriterion3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "TurbulenceFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {
		  
		  class LESProcessingData;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This strategy class computes the Qcriterion in 3D
 * and stores it in the socket named "Qcriterion"
 * 
 * @author Willem Deconinck
 */
class Qcriterion3D : public TurbulenceFunction {

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
  Qcriterion3D(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~Qcriterion3D();

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
    return "Qcriterion3D";
  }
  
  virtual void setSocketName()
  {
    m_socketName = "Qcriterion";
  }

  virtual void compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell);

protected:
  
  /// Compute the magnitude of the vorticity vector
  CFreal getAbsoluteVorticity(const std::vector<RealVector>& gradients) const;
  
  /// Compute the magnitude of the strain rate tensor
  CFreal getAbsoluteStrainRate(const std::vector<RealVector>& gradients) const;

}; // end of class Qcriterion3D

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_Qcriterion3D_hh
