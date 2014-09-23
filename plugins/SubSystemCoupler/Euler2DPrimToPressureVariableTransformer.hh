#ifndef COOLFluiD_Numerics_SubSystemCoupler_Euler2DPrimToPressureVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Euler2DPrimToPressureVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler2DPrimToPressure transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DPrimToPressureVariableTransformer : public PreVariableTransformer {
public:

  /**
   * Default Primtructor without arguments
   */
  Euler2DPrimToPressureVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DPrimToPressureVariableTransformer();

  /**
   * Transform a state into another one
   */
  RealVector* preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original);

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
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Return the size of the transformed vector
   */
  CFuint getTransformedSize(const CFuint size)
  {
    return 1;
  }

}; // end of class Euler2DPrimToPressureVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Euler2DPrimToPressureVariableTransformer_hh
