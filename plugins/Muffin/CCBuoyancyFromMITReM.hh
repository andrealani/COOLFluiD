#ifndef COOLFluiD_Muffin_CCBuoyancyFromMITReM_hh
#define COOLFluiD_Muffin_CCBuoyancyFromMITReM_hh

#include "Muffin/CC.hh"

namespace COOLFluiD {
  namespace Muffin {

/// Coupling condition, Navier-Stokes buoyancy from species concentration
class CCBuoyancyFromMITReM : public CC {

 public:  // functions

  /// Coupling condition constructor
  CCBuoyancyFromMITReM(const std::string& name);

  /// Coupling condition destructor
  ~CCBuoyancyFromMITReM() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Necessary instructions before destructor call
  void unsetup();

  /// Apply coupling condition
  void apply(const Common::SafePtr< System > applyto);


#if 0
 protected:  // sockets functions, sinks and sources

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_rhs);
    r.push_back(&s_nodes);
    r.push_back(&s_states);
    r.push_back(&s_faceneighcell);
    return r;
  }

  /// Socket to access RHS
  Framework::DataSocketSink< CFreal > s_rhs;

  /// Socket to access nodes
  Framework::DataSocketSink< Framework::Node*, Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL > s_states;

  /// Socket to access mapping surface face to corresponding boundary element
  Framework::DataSocketSink< std::pair< CFuint,CFuint > > s_faceneighcell;
#endif


#if 0
 public:  // data (user configurable)

  // System to apply this coupling from
  std::string m_system_from;

  // System to apply this coupling to
  std::string m_system_to;
#endif

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_CCBuoyancyFromMITReM_hh

