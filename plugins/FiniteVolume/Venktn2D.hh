#ifndef COOLFluiD_Numerics_FiniteVolume_Venktn2D_hh
#define COOLFluiD_Numerics_FiniteVolume_Venktn2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Limiter.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements Venkatakrishnan limiter in 2D for FVM
 *
 * @author Mehmet Sarp Yalim
 *
 */
class Venktn2D : public Framework::Limiter<CellCenterFVMData> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  Venktn2D(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~Venktn2D();
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

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
    result.push_back(&socket_uY);
    return result;
  }

  /**
   * Compute the limiter in the current face
   */
  virtual void limit(const std::vector<std::vector<Framework::Node*> >& coord,
		     Framework::GeometricEntity* const cell,
		     CFreal* limiterValue);

  /**
   * Set up the private data
   */
  virtual void setup();
  
protected:
  
  /**
   * Compute the denominator of the limiter argument
   */
  void computeDeltaMin(const Framework::Node& coord, const Framework::State& state, CFuint iVar)
  {
    const CFuint stateID = state.getLocalID();
    const RealVector& stateCoord = state.getCoordinates();
    _deltaMin = (socket_uX.getDataHandle())(stateID,iVar,state.size())*(coord[XX] - stateCoord[XX]) + 
      (socket_uY.getDataHandle())(stateID,iVar,state.size())*(coord[YY] - stateCoord[YY]);
  }
  
protected:

  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// denominator of limiter argument 
  CFreal _deltaMin;
  
  /// user defined state vector with order of magnitude of the solution 
  std::vector<CFreal> _magnitudeValues;
  
  /// user defined coefficient for the epsilon
  CFreal _coeffEps;

  /// user defined characteristic solution length in the smooth flow region
  CFreal _length;
  
}; // end of class Venktn2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Venktn2D_hh
