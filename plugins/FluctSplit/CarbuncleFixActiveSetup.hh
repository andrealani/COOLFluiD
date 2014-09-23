#ifndef COOLFluiD_Numerics_FluctSplit_CarbuncleFixActiveSetup_hh
#define COOLFluiD_Numerics_FluctSplit_CarbuncleFixActiveSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to create the socket that will hold
/// the flag telling if the carbuncle fix is active on a
/// cell.
/// To be used in conjunction with CarbuncleFixActiveUnSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino
class FluctSplit_API CarbuncleFixActiveSetup : public FluctuationSplitCom {
public: // member functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit CarbuncleFixActiveSetup(const std::string& name);

  /// Destructor
  virtual ~CarbuncleFixActiveSetup();

  /// Configure the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Executes this command
  virtual void execute();

private: // data

  /// socket with theta clending coefficient
  Framework::DataSocketSource<CFreal> socket_fix_active;

  /// Maximum number of sub-elements for high-order computations
  CFuint m_maxsubelems;

}; // class CarbuncleFixActiveSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CarbuncleFixActiveSetup_hh

