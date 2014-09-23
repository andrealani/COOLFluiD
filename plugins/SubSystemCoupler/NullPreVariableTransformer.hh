#ifndef COOLFluiD_Numerics_SubSystemCoupler_NullPreVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_NullPreVariableTransformer_hh

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
class NullPreVariableTransformer : public PreVariableTransformer {
public:

  /**
   * Default constructor without arguments
   */
  NullPreVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  ~NullPreVariableTransformer();

  /**
   * Transform a state into another one
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original)
  {
    _transVector.resize(original.size());
    _transVector = original;
    return (&_transVector);
  }

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
    return size;
  }

}; // end of class NullPreVariableTransformer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_NullPreVariableTransformer_hh
