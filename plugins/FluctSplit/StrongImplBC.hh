#ifndef COOLFluiD_Numerics_FluctSplit_StrongImplBC_hh
#define COOLFluiD_Numerics_FluctSplit_StrongImplBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong (Dirichlet-like) implicit BC
/// @author Andrea Lani
class FluctSplit_API StrongImplBC : public FluctSplit::FluctuationSplitCom {

public:

  /// Constructor.
  StrongImplBC(const std::string& name);

  /// Default destructor
  virtual ~StrongImplBC();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Set up private data
  virtual void setup();

protected:

  /// Execute on a set of dofs
  virtual void executeOnTrs() = 0;

protected:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of arrays of flags specifying if the
  // time jacobian contribution of certain variables in boundary states
  // have to be discarded
  Framework::DataSocketSink<std::vector<bool> > socket_discardTimeJacob;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<std::valarray<Framework::State*> >
  socket_bStatesNeighbors;

  /// flag telling if initialization is performed
  bool _useForInitialization;

}; // end of class StrongImplBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongImplBC_hh
