#ifndef COOLFluiD_Muffin_ComUpdateGasRate_hh
#define COOLFluiD_Muffin_ComUpdateGasRate_hh

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/Node.hh"
#include "Framework/State.hh"
#include "Muffin/MuffinData.hh"
#include "Muffin/SystemMITReM.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Command for updating gas production rate
class ComUpdateGasRate : public MuffinCom {

 public:  // core functions

  /// Gas production rate update command constructor
  ComUpdateGasRate(const std::string& name);

  /// Gas production rate update command destructor
  ~ComUpdateGasRate() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gas production rate update (virtual implementation)
  void execute();


 private:  // sockets

  /// Socket to access nodes
  Framework::DataSocketSink < Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access surface gas tracking (per boundary TRS, per boundary element)
  Framework::DataSocketSink< std::vector< GasOnSurface > > s_gasonsurf;


 public:  // sockets

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_nodes);
    r.push_back(&s_gasonsurf);
    return r;
  }


 private:  // data

  /// MITReM system to calculate gas production rate, name (configurable)
  std::string m_mitrem_str;

  /// MITReM system to calculate gas production rate, pointer
  Common::SafePtr< SystemMITReM > m_mitrem;

  /// Vector of gas-producing reactions labels (configurable)
  std::vector< std::string > m_greactions_str;

  /// Vector of gas-producing reactions indices
  std::vector< CFuint > m_greactions;

  /// If this command has been setup
  bool m_issetup;

  /// If this command saves a file with the gas production rate
  bool m_issave;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_ComUpdateGasRate_hh
