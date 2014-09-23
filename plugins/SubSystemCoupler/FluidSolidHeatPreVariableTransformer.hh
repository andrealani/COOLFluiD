#ifndef COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPreVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPreVariableTransformer_hh

//////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/DataSocketSink.hh"
#include "Common/Trio.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
      class NavierStokesVarSet;
    }
  }

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a FluidSolidHeatPre transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */

class FluidSolidHeatPreVariableTransformer : public PreVariableTransformer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  FluidSolidHeatPreVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~FluidSolidHeatPreVariableTransformer();

  /**
   * Configuration
   */
  void configure ( Config::ConfigArgs& args )
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
    cf_assert(size == 4);
    return 2;
  }

private:

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVar;

  /// arrray of gradients
  std::vector<RealVector*> _gradients;

  //normal to the face
  RealVector _normal;

  //coordinate
  std::vector<RealVector> _coords;

  ///Value for h
  CFreal _hConst;

}; // end of class FluidSolidHeatPreVariableTransformer

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPreVariableTransformer_hh
