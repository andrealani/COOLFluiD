#ifndef COOLFluiD_Numerics_SubSystemCoupler_IncompEuler2DdPuvtToTempTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_IncompEuler2DdPuvtToTempTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Incompressible Euler2D to Temperature transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class IncompEuler2DdPuvtToTempTransformer : public PreVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  IncompEuler2DdPuvtToTempTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~IncompEuler2DdPuvtToTempTransformer();

  /**
   * Configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * set up
   */
  void setup();

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
    return 1;
  }

private:

  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _updateVar;

}; // end of class IncompEuler2DdPuvtToTempTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_IncompEuler2DdPuvtToTempTransformer_hh
