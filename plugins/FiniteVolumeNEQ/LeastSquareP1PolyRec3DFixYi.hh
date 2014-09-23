#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec3DFixYi_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec3DFixYi_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/LeastSquareP1PolyRec3D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a least square polynomial reconstructor in 3D for FVM
 *
 * @author Andrea Lani
 *
 */
class LeastSquareP1PolyRec3DFixYi : public LeastSquareP1PolyRec3D {
public:

  /**
   * Constructor
   */
  LeastSquareP1PolyRec3DFixYi(const std::string& name);

  /**
   * Default destructor
   */
  ~LeastSquareP1PolyRec3DFixYi();
  
  /**
   * Set up the private data
   */
  void setup();
  
protected:
  
  /**
   * Compute the flux in the current face
   */
  void extrapolateImpl(Framework::GeometricEntity* const face);
  
  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face, CFuint iVar,
		       CFuint leftOrRight);
  
private:
  
  /// mass fractions
  RealVector m_yiL;
  
  /// mass fractions
  RealVector m_yiR;
  
}; // end of class LeastSquareP1PolyRec3DFixYi

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec3DFixYi_hh
