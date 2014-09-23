#ifndef COOLFluiD_Muffin_CCMomentumFromBubbles_hh
#define COOLFluiD_Muffin_CCMomentumFromBubbles_hh

#include "Muffin/CC.hh"

namespace COOLFluiD {
  namespace Muffin {

/// Coupling condition, Navier-Stokes external momentum from dispersed phase
class CCMomentumFromBubbles : public CC {

 public:  // functions

  /// Coupling condition constructor
  CCMomentumFromBubbles(const std::string& name);

  /// Coupling condition destructor
  ~CCMomentumFromBubbles() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Apply coupling condition
  void apply(const Common::SafePtr< System > applyto);


 protected:  // sockets functions, sinks and sources

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(CC::needsSockets());
    r.push_back(&s_mn_volume);
    return r;
  }

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_mn_volume;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_CCMomentumFromBubbles_hh

