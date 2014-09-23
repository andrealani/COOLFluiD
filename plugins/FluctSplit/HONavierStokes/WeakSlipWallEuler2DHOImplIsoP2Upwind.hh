#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DHOImplIsoP2Upwind_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DHOImplIsoP2Upwind_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/HONavierStokes/WeakSlipWall2DHOImplIsoP2Upwind.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an implicit weak slip wall bc (E. V. Weide) for Euler2D
 *
 * @author Andrea Lani
 *
 */
class WeakSlipWallEuler2DHOImplIsoP2Upwind : public WeakSlipWall2DHOImplIsoP2Upwind {
public:

  /**
   * Constructor.
   */
  WeakSlipWallEuler2DHOImplIsoP2Upwind(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEuler2DHOImplIsoP2Upwind();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:

  /**
   * Compute the normal flux and the corresponding jacobian
   */
  void computeNormalFluxAndJacob(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux,
				 RealMatrix& fluxJacob);
protected:

  /// Euler var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;


}; // end of class WeakSlipWallEuler2DHOImplIsoP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DHOImplIsoP2_hh
