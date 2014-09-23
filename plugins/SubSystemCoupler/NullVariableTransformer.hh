#ifndef COOLFluiD_Numerics_SubSystemCoupler_NullVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_NullVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "PostVariableTransformer.hh"

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
class NullVariableTransformer : public PostVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  NullVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~NullVariableTransformer();

  /**
   * Transform a vector into another one
   */
  virtual RealVector* transform(const std::vector<GeoEntityIdx>& faces,
                                const RealVector& coord,
                                const RealVector& original,
                                const RealVector& pastTransformedVector)
  {
    _transVector.resize(original.size());
    _transVector = original;

    return (&_transVector);
  }

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

}; // end of class NullVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_NullVariableTransformer_hh
