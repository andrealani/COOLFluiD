#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD3DProjectionImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD3DProjectionImpl_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak slip wall bc for MHD3DProjection
 *
 * @author Radka Keslerova
 *
 */

class WeakSlipWallMHD3DProjectionImpl : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallMHD3DProjectionImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallMHD3DProjectionImpl();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute the normal flux and the corresponding jacobian
   */
  void computeNormalFluxAndJacob(const Framework::State& state,
                                 const RealVector& normal,
                                 RealVector& flux,
                                 RealMatrix& fluxJacob);

  /**
   * Compute the face normal
   */
  void setFaceNormal(const CFuint faceID,
                     RealVector& normal);

protected:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// MHD var set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// physical data
  RealVector _physicalData;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im0;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in0;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im1;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in1;

  /// says if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;

  /// temporary flux jacobian matrix
  RealMatrix _tJacob;

  /// phi value that is to be fixed
  CFreal _refPhi;

  /// distribution coefficient
  CFreal _alpha;

}; // end of class WeakSlipWallMHD3DProjectionImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallMHD3DProjectionImpl_hh
