#ifndef COOLFluiD_Numerics_FluctSplitNEQ_ExtraDissActiveSetup_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_ExtraDissActiveSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to create the socket that will hold
/// the flag telling if the carbuncle fix is active on a
/// cell.
/// To be used in conjunction with ExtraDissActiveUnSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino
class ExtraDissActiveSetup : public FluctSplit::FluctuationSplitCom {//FluctSplitNEQ_API 
public: // member functions
  
  /// Constructor
  explicit ExtraDissActiveSetup(const std::string& name);

  /// Destructor
  virtual ~ExtraDissActiveSetup();

  /// Configure the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Executes this command
  virtual void execute();

private: // data

  /// socket with theta clending coefficient
  Framework::DataSocketSource<CFreal> socket_ExtraDiss_active;

}; // class ExtraDissActiveSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_ExtraDissActiveSetup_hh

