#ifndef COOLFluiD_Numerics_FluctSplit_ComputeVariablesDerivatives_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeVariablesDerivatives_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class computes the derivative of the variable
/// @author Thomas Wuilbaut
class FluctSplit_API ComputeVariablesDerivatives : public Framework::DataProcessingCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ComputeVariablesDerivatives(const std::string& name);

  /// Default destructor
  ~ComputeVariablesDerivatives();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Execute on a set of dofs
  void execute();

  /// Configures this object with supplied arguments.
  void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


private: //data

  /// socket for the derivatives storage
  Framework::DataSocketSource<RealVector> socket_variablesDerivatives;

  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage for flag
  Framework::DataSocketSource<CFuint> socket_flagStates;

}; // end of class ComputeVariablesDerivatives

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeVariablesDerivatives_hh
