#ifndef COOLFluiD_Numerics_LESDataProcessing_CFL_hh
#define COOLFluiD_Numerics_LESDataProcessing_CFL_hh

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
 * and stores it in a socket named "CFL".
 * 
 * @author Willem Deconinck
 */
class CFL : public TurbulenceFunction {

public: // functions
    
	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  CFL(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~CFL();

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
    return "CFL";
  }
  
  virtual void setSocketName()
  {
    m_socketName = "CFL";
  }

  virtual void compute(const RealVector& states, std::vector<RealVector>& gradients, const CFuint& iCell);

protected:
  
  /// Timestep to base CFL on
  CFreal m_dt;
  
  /// updateCoeff datahandle
  Framework::DataHandle<CFreal> m_updateCoeff;
  
  /// dimensional state
  RealVector m_primState;
  
  
}; // end of class CFL

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_CFL_hh
