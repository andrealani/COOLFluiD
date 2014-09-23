#ifndef COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPostVariableTransformerFVMCC_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPostVariableTransformerFVMCC_hh

//////////////////////////////////////////////////////////////////////

#include "PostVariableTransformer.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
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

class FluidSolidHeatPostVariableTransformerFVMCC : public PostVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  FluidSolidHeatPostVariableTransformerFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~FluidSolidHeatPostVariableTransformerFVMCC();

  /**
   * Configuration
   */
  void configure ( Config::ConfigArgs& args )
  {
    PostVariableTransformer::configure(args);
  }

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
    cf_assert(size == 1);
    return 1;
  }

private:

  /// corresponding diffusive variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSet;


}; // end of class FluidSolidHeatPostVariableTransformerFVMCC

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPostVariableTransformerFVMCC_hh
