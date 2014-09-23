#ifndef COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToIncompEuler2DdPuvtTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToIncompEuler2DdPuvtTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PostVariableTransformer.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Incompressible to Compressible Euler2D transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DConsToIncompEuler2DdPuvtTransformer : public PostVariableTransformer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  Euler2DConsToIncompEuler2DdPuvtTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DConsToIncompEuler2DdPuvtTransformer();

  /**
   * Configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * set up
   */
  void setup();

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
    return 4;
  }

private:

  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _updateVar;

  /// thermodynamic pressure infinity
  CFreal _p0Inf;

  /// dimensional R
  CFreal _RDim;

  /// specific heat ratio
  CFreal _gamma;

}; // end of class Euler2DConsToIncompEuler2DdPuvtTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToIncompEuler2DdPuvtTransformer_hh
