#ifndef COOLFluiD_Numerics_FluctSplit_CarbuncleFixActiveUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_CarbuncleFixActiveUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "Framework/ProxyDofIterator.hh"
#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to deallocate the socket telling if 
///  the carbuncle fix is active on a given cell.
/// To be used in conjunction with CarbuncleFixActiveSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino

class FluctSplit_API CarbuncleFixActiveUnSetup : public FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit CarbuncleFixActiveUnSetup(const std::string& name);

  /// Destructor.
  virtual ~CarbuncleFixActiveUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /// Executes this command
  virtual void execute();

protected: // data

  /// socket for theta blending coefficient
  Framework::DataSocketSink<CFreal> socket_fix_active;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh

