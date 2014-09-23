#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DHOP3Impl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DHOP3Impl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/HONavierStokes/WeakSlipWall2DHOImpl.hh"

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
 * This is the discretization for P3 elements
 * It is not working.....do not know the reason why
 *
 * @author Nadege Villedieu
 *
 */
class WeakSlipWallEuler2DHOP3Impl : public WeakSlipWall2DHOImpl {
public:

  /**
   * Constructor.
   */
  WeakSlipWallEuler2DHOP3Impl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEuler2DHOP3Impl();

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
  
  
}; // end of class WeakSlipWallEuler2DHOP3Impl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DHOP3Impl_hh
