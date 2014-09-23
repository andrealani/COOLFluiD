#ifndef COOLFluiD_Numerics_LESDataProcessing_LocalTimeStep_hh
#define COOLFluiD_Numerics_LESDataProcessing_LocalTimeStep_hh

//////////////////////////////////////////////////////////////////////////////

#include "TurbulenceFunction.hh"
#include "LES/LESVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {
		  
		  class LESProcessingData;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This strategy class computes the DYNAMIC Subgrid-scale viscosity
 * and stores it in a socket named "LocalTimeStep".
 * 
 * @author Willem Deconinck
 */
class LocalTimeStep : public TurbulenceFunction {

public: // functions
    
	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  LocalTimeStep(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~LocalTimeStep();

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
    return "LocalTimeStep";
  }
  
  virtual void setSocketName()
  {
    m_socketName = "LocalTimeStep";
  }

  virtual void compute(const RealVector& states, std::vector<RealVector>& gradients, const CFuint& iCell);

protected:
    
  /// updateCoeff datahandle
  Framework::DataHandle<CFreal> m_updateCoeff;
  
  /// dimensional state
  RealVector m_primState;
  
  
}; // end of class LocalTimeStep

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_LocalTimeStep_hh
