#ifndef COOLFluiD_Numerics_LESDataProcessing_Qcriterion2D_hh
#define COOLFluiD_Numerics_LESDataProcessing_Qcriterion2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "TurbulenceFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {
		  
		  class LESProcessingData;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This strategy class computes the Qcriterion in 2D
 * and stores it in the socket named "Qcriterion"
 * 
 * @author Willem Deconinck
 */
class Qcriterion2D : public TurbulenceFunction {

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
  Qcriterion2D(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~Qcriterion2D();

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
    return "Qcriterion2D";
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

}; // end of class Qcriterion2D

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_Qcriterion2D_hh
