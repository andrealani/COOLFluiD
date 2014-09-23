#ifndef COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNS2DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNS2DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler2D
/// @author Andrea Lani
class FluctSplit_API StrongNoSlipWallAdiabaticNS2DImpl : public FluctuationSplitCom {

public:

  /// Constructor.
  StrongNoSlipWallAdiabaticNS2DImpl(const std::string& name);

  /// Default destructor
  ~StrongNoSlipWallAdiabaticNS2DImpl();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Execute on a set of dofs
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

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<std::valarray<Framework::State*> > socket_bStatesNeighbors;

  /// flag telling if initialization is performed
  bool _useForInitialization;

}; // end of class StrongNoSlipWallAdiabaticNS2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNS2DImpl_hh
