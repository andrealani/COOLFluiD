#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DTurb_hh

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
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 */
class LeastSquareP1PolyRec2DTurb : public LeastSquareP1PolyRec2D {
public:

  /**
   * Constructor
   */
  LeastSquareP1PolyRec2DTurb(const std::string& name);

  /**
   * Default destructor
   */
  ~LeastSquareP1PolyRec2DTurb();
  
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

}; // end of class LeastSquareP1PolyRec2DTurb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DTurb_hh
