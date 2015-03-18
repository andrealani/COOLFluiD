#ifndef COOLFluiD_Numerics_FiniteVolume_SA2DDES97_hh
#define COOLFluiD_Numerics_FiniteVolume_SA2DDES97_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeTurb/SA2DSourceTerm.hh" 

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the implementation of the SA - DES97 model proposed by the P. R. Spalart for 2D.
 * @see paper "COMMENTS OF THE FEASIBILITY OF LES FOR WINGS, AND ON A HYBRID RANS/LES APPROACH",
 * Proceedings of the first AFOSR International Conference on DNS/LES".
 * 
 * @author Christos Gkoudesnes
 *
 */
class SA2DDES97 : public SA2DSourceTerm {

public:

  /**
   * Constructor
   * @see SA2DSourceTerm
   */
  SA2DDES97(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SA2DDES97();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    SA2DSourceTerm ::configure(args);
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
    SA2DSourceTerm ::computeSource( element, source, jacobian);
  }
  
protected:  
  
  // method that should override the RANS method
  ///@return the DES97 distance 
  virtual CFreal getDistance (Framework::GeometricEntity *const element);
  
  ///@return the Subgrid length-scale for the DES modes multiplying by the DES constant
  virtual CFreal getSubLenScale (Framework::GeometricEntity *const element);
  
}; // end of class SA2DDES97

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SA2DDES97_hh
