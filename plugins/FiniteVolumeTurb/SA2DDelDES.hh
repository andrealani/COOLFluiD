#ifndef COOLFluiD_Numerics_FiniteVolume_SA2DDelDES_hh
#define COOLFluiD_Numerics_FiniteVolume_SA2DDelDES_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeTurb/SA2DDES97.hh" 

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the implementation of the SA - DDES model proposed by the P. R. Spalart et all for 2D.
 * @see paper "A New Version of DES, Resistant to Ambiguouis Grid Densities", Theoretical and Computational Fluid Dynamics,
 * Vol. 20, 2006, pp. 181-195.
 * 
 * @author Christos Gkoudesnes
 *
 */
class SA2DDelDES : public SA2DDES97 {

public:

  /**
   * Constructor
   * @see SA2DDES97
   */
  SA2DDelDES(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SA2DDelDES();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    SA2DDES97::configure(args);
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
   /**
   * Compute the source term
   */
 virtual void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian)
 
  {
    SA2DDES97 ::computeSource( element, source, jacobian);
  }
   
protected:  
  
  // method that should override the DES97 method
  ///@return the DDES distance to the wall
  virtual CFreal getDistance (Framework::GeometricEntity *const element);
  
  // method that compute the delaying function (fd) of the DDES mode
  ///@return the fd factor
  virtual CFreal computeDelayingFunction (Framework::GeometricEntity *const element);
  
}; // end of class SA2DDelDES

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SA2DDelDES_hh
