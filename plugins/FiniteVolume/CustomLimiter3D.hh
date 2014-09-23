#ifndef COOLFluiD_Numerics_FiniteVolume_CustomLimiter3D_hh
#define COOLFluiD_Numerics_FiniteVolume_CustomLimiter3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "CustomLimiter2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a family of user defined limiters associated to 3D
 * high-order reconstruction
 *
 * @author Andrea Lani
 *
 */
class CustomLimiter3D : public CustomLimiter2D {
public:
  
  /**
   * Constructor
   */
  CustomLimiter3D(const std::string& name);
  
  /**
   * Default destructor
   */
   ~CustomLimiter3D();
  
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
      CustomLimiter2D::needsSockets();
    result.push_back(&socket_uZ);
    
    return result;
  }
  
protected:
  
  /**
   * Compute the denominator of the limiter argument
   */
  void computeDeltaMin(const Framework::Node& coord,
		       const Framework::State& state, CFuint iVar);
  
protected:
  
  /// socket for uZ values
  Framework::DataSocketSink<CFreal> socket_uZ;
  
}; // end of class CustomLimiter3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CustomLimiter3D_hh
