#ifndef COOLFluiD_Numerics_LESDataProcessing_SGSViscosity_hh
#define COOLFluiD_Numerics_LESDataProcessing_SGSViscosity_hh

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
 * and stores it in a socket named "SGSViscosity".
 * 
 * @author Willem Deconinck
 */
class SGSViscosity : public TurbulenceFunction {

public: // functions
  
  typedef LES::LESVarSet LESVAR;
  
	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  SGSViscosity(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~SGSViscosity();

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
    return "SGSViscosity";
  }
  
  virtual void setSocketName()
  {
    m_socketName = "SGSViscosity";
  }

  virtual void compute(const RealVector& states, std::vector<RealVector>& gradients, const CFuint& iCell);

  virtual void computeAverage(const CFreal& oldWeight, const CFreal& newWeight, const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell);
protected:
  
  /// diffusive varSet
  Common::SafePtr<LESVAR> m_lesVar;
  
  /// number of equations
  CFuint m_nbEqs;

  /// Conversion interface for gradients
  std::vector<RealVector*> m_gradientPtrs;
  
  /// dimensional state
  RealVector m_dimState;
  
}; // end of class SGSViscosity

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_SGSViscosity_hh
