#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DBcFix_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DBcFix_hh

//////////////////////////////////////////////////////////////////////////////

#include "LeastSquareP1PolyRec2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a least square polynomial reconstructor in 2D for FVM
 *
 * @author Andrea Lani
 */
class LeastSquareP1PolyRec2DBcFix : public LeastSquareP1PolyRec2D {
public:

  /**
   * Constructor
   */
  LeastSquareP1PolyRec2DBcFix(const std::string& name);

  /**
   * Default destructor
   */
  ~LeastSquareP1PolyRec2DBcFix();

protected:

  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face);
  
  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face,
		       CFuint iVar, CFuint leftOrRight);
private:
  
  // distance between inner and ghost state
  CFreal m_drXiXg;
  
  // distance between quadrature point and ghost state
  CFreal m_drXqXg;

  // distance between quadrature point and inner state
  CFreal m_drXqXi;
  
}; // end of class LeastSquareP1PolyRec2DBcFix

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DBcFix_hh
