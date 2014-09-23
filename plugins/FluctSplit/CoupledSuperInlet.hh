#ifndef COOLFluiD_Numerics_FluctSplit_CoupledSuperInlet_hh
#define COOLFluiD_Numerics_FluctSplit_CoupledSuperInlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluctuationSplitData.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a supersonic inlet command for use in coupled subsystems
/// @author Thomas Wuilbaut
class FluctSplit_API CoupledSuperInlet : public FluctuationSplitCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  CoupledSuperInlet(const std::string& name);

  /// Default destructor
  ~CoupledSuperInlet();

  /// Configures the command.
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Execute on the current TRS
  void executeOnTrs();

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the isUpdated flag
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// name of the datahandle containing the nodal values
  std::string _interfaceName;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // end of class CoupledSuperInlet

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CoupledSuperInlet_hh
