#ifndef COOLFluiD_Numerics_FiniteVolume_CustomLimiter2D_hh
#define COOLFluiD_Numerics_FiniteVolume_CustomLimiter2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "CustomLimiter1D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a family of user defined limiters associated to 2D
 * high-order reconstruction
 *
 * @author Andrea Lani
 *
 */
class CustomLimiter2D : public CustomLimiter1D {
public:

  /**
   * Constructor
   */
  CustomLimiter2D(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~CustomLimiter2D();
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
      CustomLimiter1D::needsSockets();
    result.push_back(&socket_uY);
    
    return result;
  }

protected:
  
  /**
   * Compute the denominator of the limiter argument
   */
  virtual void computeDeltaMin(const Framework::Node& coord, 
			       const Framework::State& state, CFuint iVar);
  
protected:
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;
  
}; // end of class CustomLimiter2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CustomLimiter2D_hh
