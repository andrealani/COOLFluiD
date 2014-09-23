#ifndef COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressureVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressureVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"
#include "NavierStokes/EulerPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler2DConsToPressure transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DConsToPressureVariableTransformer : public PreVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DConsToPressureVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DConsToPressureVariableTransformer();

  /**
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Transform a state into another one
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
		       const RealVector& coord,
		       const RealVector& original);

  /**
   * Transform a vector into another one
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
    return 2;
  }

protected:

  ///Link to the Heat Physical Model
  Common::SafePtr<Physics::NavierStokes::EulerPhysicalModel<DIM_2D> > _model;


}; // end of class Euler2DConsToPressureVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressureVariableTransformer_hh
