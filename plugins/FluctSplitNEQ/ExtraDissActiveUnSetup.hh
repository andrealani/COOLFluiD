#ifndef COOLFluiD_Numerics_FluctSplitNEQ_ExtraDissActiveUnSetup_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_ExtraDissActiveUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "Framework/ProxyDofIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to deallocate the socket telling if 
///  the carbuncle fix is active on a given cell.
/// To be used in conjunction with ExtraDissActiveSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino

class ExtraDissActiveUnSetup : public FluctSplit::FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit ExtraDissActiveUnSetup(const std::string& name);

  /// Destructor.
  virtual ~ExtraDissActiveUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /// Executes this command
  virtual void execute();

protected: // data

  /// socket for theta blending coefficient
  Framework::DataSocketSink<CFreal> socket_ExtraDiss_active;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_ExtraDissActiveUnSetup_hh

