#ifndef COOLFluiD_Numerics_FiniteVolume_CustomLimiter1D_hh
#define COOLFluiD_Numerics_FiniteVolume_CustomLimiter1D_hh

//////////////////////////////////////////////////////////////////////////////

#include "CustomLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a family of user defined limiters associated to 1D
 * high-order reconstruction
 *
 * @author Andrea Lani
 *
 */
class CustomLimiter1D : public CustomLimiter {
public:

  /**
   * Constructor
   */
  CustomLimiter1D(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~CustomLimiter1D();
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
      Framework::Limiter<CellCenterFVMData>::needsSockets();
    result.push_back(&socket_stencil);
    result.push_back(&socket_uX);
    return result;
  }

  /**
   * Compute the limiter in the current face
   */
  virtual void limit(const std::vector<std::vector<Framework::Node*> >& coord,
		     Framework::GeometricEntity* const cell,
		     CFreal* limiterValue);
  
protected:
  
  /**
   * Compute the denominator of the limiter argument
   */
  virtual void computeDeltaMin(const Framework::Node& coord, 
			       const Framework::State& state, CFuint iVar);
  
protected:
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
      
}; // end of class CustomLimiter1D

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace FiniteVolume

  } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CustomLimiter1D_hh
