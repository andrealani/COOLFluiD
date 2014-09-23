#ifndef COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNSKOmega2D_hh
#define COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNSKOmega2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler2D
/// @author Andrea Lani
class FluctSplit_API StrongNoSlipWallAdiabaticNSKOmega2D : public FluctuationSplitCom {

public: // functions

  /// Constructor.
  StrongNoSlipWallAdiabaticNSKOmega2D(const std::string& name);

  /// Default destructor
  ~StrongNoSlipWallAdiabaticNSKOmega2D();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

protected: // functions

  /// Execute on a set of dofs
  void executeOnTrs();

private: // data

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// x-coponent of the wall velocity
  CFreal m_u0;

  /// y-coponent of the wall velocity
  CFreal m_v0;



}; // end of class StrongNoSlipWallAdiabaticNSKOmega2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNSKOmega2D_hh
