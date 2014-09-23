#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD2DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD2DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakSlipWall2DImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD2DVarSet;
    }
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an implicit weak slip wall bc for MHD2D
 *
 * @author Radka Keslerova
 *
 */
class WeakSlipWallMHD2DImpl : public WeakSlipWall2DImpl {
public:

  /**
   * Constructor.
   */
  WeakSlipWallMHD2DImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallMHD2DImpl();

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

  /// MHD var set
  Common::SafePtr<Physics::MHD::MHD2DVarSet> _varSet;
  
}; // end of class WeakSlipWallMHD2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD2DImpl_hh
