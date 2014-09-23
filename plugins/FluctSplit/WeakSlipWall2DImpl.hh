#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWall2DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWall2DImpl_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an implicit weak slip wall bc (E. V. Weide)
/// @author Andrea Lani
class FluctSplit_API WeakSlipWall2DImpl : public FluctuationSplitCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  WeakSlipWall2DImpl(const std::string& name);

  /// Default destructor
  virtual ~WeakSlipWall2DImpl();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Execute on a set of dofs
  virtual void executeOnTrs();

  /// Compute the normal flux and the corresponding jacobian
  virtual void computeNormalFluxAndJacob(const Framework::State& state,
  				 const RealVector& normal,
  				 RealVector& flux,
  				 RealMatrix& fluxJacob) = 0;

  /// Compute the face normal
  void setFaceNormal(const CFuint faceID,
      RealVector& normal);

protected:

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

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

  /// distribution coefficient
  CFreal m_alpha;

  /// vector of fluxes at boundary
  std::vector<RealVector> m_flux;

  /// vector of Jacobian fluxes at boundary
  std::vector<RealMatrix> m_fluxJacob;

  /// vector of states
  std::vector<Framework::State*> m_state;


}; // end of class WeakSlipWall2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWall2DImpl_hh
