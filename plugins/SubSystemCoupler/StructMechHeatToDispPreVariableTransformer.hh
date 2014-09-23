#ifndef COOLFluiD_Numerics_SubSystemCoupler_StructMechHeatToDispPreVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StructMechHeatToDispPreVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a null transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMechHeatToDispPreVariableTransformer : public PreVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  StructMechHeatToDispPreVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMechHeatToDispPreVariableTransformer();

  /**
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Transform a state into another one
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original);

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

}; // end of class StructMechHeatToDispPreVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StructMechHeatToDispPreVariableTransformer_hh
