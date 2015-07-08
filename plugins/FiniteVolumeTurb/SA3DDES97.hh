#ifndef COOLFluiD_Numerics_FiniteVolume_SA3DDES97_hh
#define COOLFluiD_Numerics_FiniteVolume_SA3DDES97_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeTurb/SA3DSourceTerm.hh" 
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the implementation of the SA - DES97 model proposed by the P. R. Spalart for 3D.
 * @see paper "COMMENTS OF THE FEASIBILITY OF LES FOR WINGS, AND ON A HYBRID RANS/LES APPROACH",
 * Proceedings of the first AFOSR International Conference on DNS/LES".
 * 
 * @author Christos Gkoudesnes
 *
 */
class SA3DDES97 : public SA3DSourceTerm {

public:

  /**
   * Constructor
   * @see SA3DSourceTerm
   */
  SA3DDES97(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SA3DDES97();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    SA3DSourceTerm ::configure(args);
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
    SA3DSourceTerm ::computeSource( element, source, jacobian);
  }
  
   /**
      Returns the DataSocket's that this command provides as sources
      @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
protected:  
  
  // method that should override the RANS method
  ///@return the DES97 distance 
  virtual CFreal getDistance (Framework::GeometricEntity *const element);
  
  ///@return the Subgrid length-scale for the DES modes multiplying by the DES constant
  virtual CFreal getSubLenScale (Framework::GeometricEntity *const element);
  

protected: // data

  /// socket for Subgrid length scale
  Framework::DataSocketSource<CFreal> socket_length_scale;  
  
  /// socket for switch function
  Framework::DataSocketSource<CFreal> socket_switch_function; 

}; // end of class SA3DDES97

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SA3DDES97_hh
