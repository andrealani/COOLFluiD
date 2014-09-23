#ifndef COOLFluiD_Numerics_SubSystemCoupler_MeshMovementPredictorVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_MeshMovementPredictorVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class MeshMovementPredictorVariableTransformer : public PreVariableTransformer {
public:

 /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  MeshMovementPredictorVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~MeshMovementPredictorVariableTransformer();

  /**
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Transform a state into another one
   */
  RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
		           const RealVector& coord,
		           const RealVector& original)
  {
    cf_assert(false);
    return 0;
  }

  /**
   * Transform a state into another one
   */
  RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
		           const RealVector& coord,
		           const RealVector& original,
                           const RealVector& shapeFunctions);

  /**
   * Returns the DataSocket's that this strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Return the size of the transformed vector
   */
  CFuint getTransformedSize(const CFuint size)
  {
    return 2;
  }

protected:

  /// handle to the pastStates, pastVel and pastAccel
  Framework::DataSocketSink<Framework::State* >  socket_pastStates;
  Framework::DataSocketSink<Framework::State* >  socket_pastStatesD;
  Framework::DataSocketSink<Framework::State* >  socket_pastStatesD2;

  /// values of the pastStates, pastVel and pastAccel
  RealVector _pastDisp;
  RealVector _pastVel;
  RealVector _pastAcc;

  CFreal _alpha0;
  CFreal _alpha1;

}; // end of class MeshMovementPredictorVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_MeshMovementPredictorVariableTransformer_hh
