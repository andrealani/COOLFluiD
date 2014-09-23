#ifndef COOLFluiD_Numerics_FluctSplit_UpdateCoeffSetup_hh
#define COOLFluiD_Numerics_FluctSplit_UpdateCoeffSetup_hh

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
/// executed in order to create the socket to the theta of blending scheme.
/// ( \Phi^B = \theta \Phi^N +(1-\theta)*\phi^LDA)
/// To be used in conjunction with ThetaUnSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino
class FluctSplit_API UpdateCoeffSetup : public FluctuationSplitCom {
public: // member functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit UpdateCoeffSetup(const std::string& name);

  /// Destructor
  virtual ~UpdateCoeffSetup();

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
  Framework::DataSocketSource<CFreal> socket_updateCoeff;

  /// Maximum number of sub-elements for high-order computations
//   CFuint m_maxsubelems;

}; // class UpdateCoeffSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UpdateCoeffSetup_hh

