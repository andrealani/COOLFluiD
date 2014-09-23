#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD2DProjectionImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD2DProjectionImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakSlipWall2DImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD2DProjectionVarSet;
    }
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak slip wall bc for MHD2DProjection
 *
 * @author Radka Keslerova
 *
 */

class WeakSlipWallMHD2DProjectionImpl : public WeakSlipWall2DImpl {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallMHD2DProjectionImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallMHD2DProjectionImpl();

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
  Common::SafePtr<Physics::MHD::MHD2DProjectionVarSet> _varSet;
  
  /// phi value that is to be fixed
  CFreal _refPhi;
  
}; // end of class WeakSlipWallMHD2DProjectionImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD2DProjectionImpl_hh
