#ifndef COOLFluiD_Muffin_CC_hh
#define COOLFluiD_Muffin_CC_hh

#include "Common/COOLFluiD.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Base class for coupling condition
class CC : public MuffinCom,
           public Common::TaggedObject {

 public:  // functions

  /// Coupling condition constructor
  CC(const std::string& name);

  /// Coupling condition destructor
  ~CC() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  virtual void setup();

  /// Apply coupling condition
  virtual void apply(const Common::SafePtr< System > applyto) = 0;


 protected:  // functions

  /// Log information, debug and error messages
  void log(const std::string& msg) { getMethodData().log("CC " + getName() + ": " + msg); }
  void ver(const std::string& msg) { getMethodData().ver("CC " + getName() + ": " + msg); }
  void err(const std::string& msg) { getMethodData().err("CC " + getName() + ": " + msg); }


 public:  // sockets functions

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_rhs);
    r.push_back(&s_nodes);
    r.push_back(&s_states);
    r.push_back(&s_faceneighcell);
    return r;
  }


 protected:  // sockets

  /// Socket to access RHS
  Framework::DataSocketSink< CFreal > s_rhs;

  /// Socket to access nodes
  Framework::DataSocketSink< Framework::Node*, Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL > s_states;

  /// Socket to access mapping surface face to corresponding boundary element
  Framework::DataSocketSink< std::pair< CFuint,CFuint > > s_faceneighcell;


 private:  // data

  // System command (name) to apply this coupling condition from (configurable)
  std::string m_applyfrom_str;


 protected:  // data

  // System command to apply this coupling condition from
  Common::SafePtr< System > m_applyfrom;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_CC_hh

