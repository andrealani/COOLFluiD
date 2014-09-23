#ifndef COOLFluiD_Numerics_FluctSplit_StrongSlipWallMHD2DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSlipWallMHD2DImpl_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for MHD2DImpl
 *
 * @author Radka Keslerova
 *
 *
 *
 */
class StrongSlipWallMHD2DImpl : public FluctuationSplitCom {

public:

  /**
   * Constructor.
   */
  StrongSlipWallMHD2DImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSlipWallMHD2DImpl();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

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

private:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the state's neighbors
  Framework::DataSocketSink<std::valarray<Framework::State*> > socket_bStatesNeighbors;

  /// socket for normals
  Framework::DataSocketSink<InwardNormalsData* > socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > _bcNormals;

  /// storage for a block of values to insert in the jacobian matrix
  RealMatrix _block;

  /// storage for a block of values got from the jacobian matrix
  RealMatrix _jacobElem;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in;

  /// storage for values to subtract from a block diagonal
  RealVector _diag;

}; // end of class StrongSlipWallMHD2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallMHD2DImpl_hh
