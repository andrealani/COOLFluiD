#ifndef COOLFluiD_Numerics_FiniteVolume_NSFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_NSFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeDiffusiveFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model
 *
 * @author Andrea Lani
 *
 */
template <typename DIFFVS>      
class NSFlux : public ComputeDiffusiveFlux {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NSFlux(const std::string& name);

  /**
   * Default destructor
   */
  ~NSFlux();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up 
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeDiffusiveFlux::needsSockets();

    result.push_back(&socket_wallDistance);
    result.push_back(&socket_volumes);
    result.push_back(&socket_faceCenters);
    
    return result;
  }
  
  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);

protected: //helper functions

  /**
   * Set the wall distance value in the diffusive var set if needed
   */
  void setWallDistance() 
  {
    using namespace COOLFluiD::Framework;
    
    const GeometricEntity& geo = *getMethodData().getCurrentFace();
    const CFreal distance = 0.5*(_wallDistance[geo.getState(0)->getLocalID()] +
				 _wallDistance[geo.getState(1)->getLocalID()]);
    _diffVar->setWallDistance(distance);
  }
  
protected: // data
  
  /// socket for the wallDistance storage
  Framework::DataSocketSink<CFreal> socket_wallDistance;
  
  /// socket for the cell volumes storage
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// socket for the face centers storage
  Framework::DataSocketSink<CFreal> socket_faceCenters;
  
  /// data handle cashed for efficiency  reasons
  Framework::DataHandle<CFreal> _wallDistance;
  
  /// update variable set
  Common::SafePtr<DIFFVS> _diffVar;
  
  /// number states in the control volume
  CFuint _nbCVStates;
  
  /// cell radius
  CFreal _radius;
  
  /// flag to tell if wall distance exists
  bool _wallDistanceExists;
  
  /// name of the wall distance data handle
  std::string _wallDistanceDataHandleName;
  
  // array of states involved in the flux computation
  std::vector<RealVector*> _states;

  // array of values (p, u, v, T, ...)
  RealMatrix _values;
  
  /// arrray of gradients
  std::vector<RealVector*> _gradients;
  
  // array of average state in update variables
  RealVector _avState;
  
  // average radius vector
  RealVector _avRadiusVec;
  
  /// flag telling if the radius is needed
  bool _isRadiusNeeded;
  
}; // end of class NSFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NSFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NSFlux_hh
