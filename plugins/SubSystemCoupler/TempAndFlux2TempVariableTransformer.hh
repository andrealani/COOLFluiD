#ifndef COOLFluiD_Numerics_SubSystemCoupler_TempAndFlux2TempVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_TempAndFlux2TempVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PostVariableTransformer.hh"
#include "Framework/GeometricEntity.hh"
#include "Heat/HeatPhysicalModel.hh"
#include "Framework/DataSocketSink.hh"
#include "Common/Trio.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a TempAndFlux2Flux transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */

class TempAndFlux2TempVariableTransformer : public PostVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  TempAndFlux2TempVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~TempAndFlux2TempVariableTransformer();

  /**
   * Configuration
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    PostVariableTransformer::configure(args);
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
   * Transform a vector into another one
   */
  virtual RealVector* transform(const std::vector<GeoEntityIdx>& faces,
                                const RealVector& coord,
                                const RealVector& original,
                                const RealVector& pastTransformedVector);

  /**
   * Transform a vector into another one (in the case of nodal values)
   */
  virtual RealVector* transform(const std::vector<GeoEntityIdx>& faces,
                                const RealVector& coord,
                                const RealVector& currentState,
                                const RealVector& original,
                                const RealVector& pastTransformedVector)
  {
    return transform(faces, coord, original, pastTransformedVector);
  }

  /**
   * Return the size of the transformed vector
   */
  CFuint getTransformedSize(const CFuint size)
  {
    cf_assert(size == 2);

    return 1;
  }

private:

  ///Link to the Heat Physical Model
  Common::SafePtr<Physics::Heat::HeatPhysicalModel> _model;


  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

}; // end of class TempAndFlux2TempVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_TempAndFlux2TempVariableTransformer_hh
