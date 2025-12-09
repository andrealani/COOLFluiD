#ifndef COOLFluiD_Numerics_FiniteVolume_Venktn3DStrict_hh
#define COOLFluiD_Numerics_FiniteVolume_Venktn3DStrict_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/Venktn2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements Venkatakrishnan limiter in 3D for FVM
 *
 * @author Mehmet Sarp Yalim
 *
 */
class Venktn3DStrict : public Venktn2D {
public:
  

  static void defineConfigOptions(Config::OptionList& options);
  /**
   * Constructor
   */
  Venktn3DStrict(const std::string& name);

  /**
   * Default destructor
   */
  ~Venktn3DStrict();

 /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
      Venktn2D::needsSockets();
    result.push_back(&socket_uZ);
    
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
  virtual void setup()
  {
    Venktn2D::setup();
  }

protected:
  
  /**
   * Compute the denominator of the limiter argument
   */
  void computeDeltaMin(const Framework::Node& coord, const Framework::State& state, CFuint iVar)
  {
    const CFuint stateID = state.getLocalID();
    const RealVector& stateCoord = state.getCoordinates();
    _deltaMin = (socket_uX.getDataHandle())(stateID,iVar,state.size())*(coord[XX] - stateCoord[XX]) + 
      (socket_uY.getDataHandle())(stateID,iVar,state.size())*(coord[YY] - stateCoord[YY]) +
      (socket_uZ.getDataHandle())(stateID,iVar,state.size())*(coord[ZZ] - stateCoord[ZZ]);
  }

protected:
  CFreal _strictCoeff;  

  bool _psiMinEqual1; 

  /// corresponding to variables which don't implement limiter
  std::vector<CFuint> _NoLimiterID;
  
private:

  /// socket for uZ values
  Framework::DataSocketSink<CFreal> socket_uZ;
  
}; // end of class Venktn3DStrict

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Venktn3DStrict_hh
