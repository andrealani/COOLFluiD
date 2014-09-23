#ifndef COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter3D_hh
#define COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter3D_hh

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
class ReconstructionFilter3D : public ReconstructionFilter {

public: // functions

  // /** 
  //  * Defines the Config Options of this class
  //  * @param    options   an OptionList where to add the options
  //  */
  //   static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  ReconstructionFilter3D(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~ReconstructionFilter3D();

  /** 
   * Set private data that will be used during the computation
   */
   virtual void setup();

  /** 
   * Configure the object
   */
  // virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ReconstructionFilter3D";
  }

  /**
   * Show the reconstruction equation
   *
   * @param   order   order of the reconstruction
   * @return  LaTeX formatted string for the equation used in the reconstruction
   */
  virtual std::string reconstructionEquation(const CFuint order) const;
  
  void outputTransferFunction();
  
protected: // helper functions
  
  /**
   * Calculate the system matrix of the least squares problem Ax=b
   * @param A the system matrix to be calculated
   * @param stencil The stencil to be used for the calculation
   */
  void calculateSystemMatrix(RealMatrix& A);
  
  /**
   * Coefficient of x^n1 y^n2 z^n3 in expansion of (x+y+z)^n (with n=n1+n2+n3)
   *
   * @param   n1   power corresponding to x-direction
   * @param   n2   power corresponding to y-direction
   * @param   n3   power corresponding to z-direction
   * @return  trinomial coefficient
   */
  CFreal trinomialCoefficient(const CFuint n1, const CFuint n2, const CFuint n3) const; 

private: // data

	  
protected: // data


}; // end of class FilterStrategy

//////////////////////////////////////////////////////////////////////////////

  	} // namespace ExplicitFilters

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter3D_hh
