#ifndef COOLFluiD_Numerics_SubSystemCoupler_TestVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_TestVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PostVariableTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Test transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class TestVariableTransformer : public PostVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  TestVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~TestVariableTransformer();

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
    return size;
  }

}; // end of class TestVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_TestVariableTransformer_hh
