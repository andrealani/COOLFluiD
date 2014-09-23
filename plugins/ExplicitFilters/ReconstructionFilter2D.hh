#ifndef COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter2D_hh
#define COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "ReconstructionFilter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

/** 
 * This class offers a basic interface for the different explicit
 * filter
 * 
 * @author Willem Deconinck
 */
class ReconstructionFilter2D : public ReconstructionFilter {

  
public: // functions

  // /** 
  //  * Defines the Config Options of this class
  //  * @param    options   an OptionList where to add the options
  //  */
  //   static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  ReconstructionFilter2D(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~ReconstructionFilter2D();

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();

  // /** 
  //  * Configure the object
  //  */
  // virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ReconstructionFilter2D";
  }
  

  /**
   * Show the reconstruction equation
   *
   * @param   order   order of the reconstruction
   * @return  LaTeX formatted string for the equation used in the reconstruction
   */
  std::string reconstructionEquation(const CFuint order) const;
  
  
protected: // helper functions

  /**
   * Coefficient of x^n1 y^n2 in expansion of (x+y)^n (with n=n1+n2)
   *
   * @param   n1   power corresponding to x-direction
   * @param   n2   power corresponding to y-direction
   * @return  trinomial coefficient
   */
  CFreal binomialCoefficient(const CFuint n1, const CFuint n2) const; 
  
  /**
   * Calculate the system matrix of the least squares problem Ax=b
   * @param A the system matrix to be calculated
   * @param stencil The stencil to be used for the calculation
   */
  void calculateSystemMatrix(RealMatrix& A);


private: // data
	  
protected: // data


}; // end of class FilterStrategy

//////////////////////////////////////////////////////////////////////////////

  	} // namespace ExplicitFilters

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter2D_hh
