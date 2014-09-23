#ifndef COOLFluiD_Numerics_SubSystemCoupler_Temp2TempAndFluxVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Temp2TempAndFluxVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"
#include "Framework/GeometricEntity.hh"
#include "Heat/Heat2D.hh"
#include "Framework/DataSocketSink.hh"
#include "Common/Trio.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Temp2TempAndFlux transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */

class Temp2TempAndFluxVariableTransformer : public PreVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Temp2TempAndFluxVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~Temp2TempAndFluxVariableTransformer();

  /**
   * Configuration
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    PreVariableTransformer::configure(args);
  }

  /**
   * Returns the DataSocket's that this strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Transform a vector
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
                                   const RealVector& coord,
                                   const RealVector& original);

  /**
   * Transform a vector
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
                                   const RealVector& coord,
                                   const RealVector& original,
                                   const RealVector& shapeFunctions)
  {
    return preTransform(faces, coord, original);
  }

  /**
   * Return the size of the transformed vector
   */
  CFuint getTransformedSize(const CFuint size)
  {
    cf_assert(size == 1);

    return 2;
  }

private:

  ///Link to the Heat Physical Model
  Common::SafePtr<Physics::Heat::HeatPhysicalModel> _model;


  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  //gradients of temperature
  RealVector _gradientsT;

  //normal to the face
  RealVector _normal;

  //coordinate
  std::vector<RealVector> _coords;

}; // end of class Temp2TempAndFluxVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Temp2TempAndFluxVariableTransformer_hh
