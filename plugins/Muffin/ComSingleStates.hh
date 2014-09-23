#ifndef COOLFluiD_Muffin_ComSingleStates_hh
#define COOLFluiD_Muffin_ComSingleStates_hh

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/Node.hh"
#include "Framework/State.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Command for queuing multiple states
class ComSingleStates : public MuffinCom {

 public:  // core functions

  /// Single states queue command constructor
  ComSingleStates(const std::string& name);

  /// Single states queue command destructor
  ~ComSingleStates() {}

  /// Defines the Config Option's of this command
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase (virtual implementation)
  void setup();

  /// Iterate on solution and derivatives queues (virtual implementation)
  void execute();


 private:  // sockets

  /// Socket to access nodes
  Framework::DataSocketSink < Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink < Framework::State*,Framework::GLOBAL > s_states;

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_nvolume;

  /// Sockets to provide solution and (nodal) spacial derivative x, y and z components
  Framework::DataSocketSource< CFreal > s_sol;
  Framework::DataSocketSource< CFreal > s_ddx;
  Framework::DataSocketSource< CFreal > s_ddy;
  Framework::DataSocketSource< CFreal > s_ddz;


 public:  // sockets

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_nodes);
    r.push_back(&s_states);
    r.push_back(&s_nvolume);
    return r;
  }

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > r;
    r.push_back(&s_sol);
    r.push_back(&s_ddx);
    r.push_back(&s_ddy);
    r.push_back(&s_ddz);
    return r;
  }


 private:  // data

  /// If solution is queued (configurable)
  CFuint m_qsolution;

  /// If solution spacial derivatives are queued (configurable)
  CFuint m_qderivatives;

  /// Number of nodes
  CFuint nbNodes;

  /// Number of states
  CFuint nbStates;

  /// Number of dimensions (Node size)
  CFuint nbDim;

  /// Number of variables per state (State size)
  CFuint nbVar;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_ComSingleStates_hh
