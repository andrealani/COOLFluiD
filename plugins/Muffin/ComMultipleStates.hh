#ifndef COOLFluiD_Muffin_ComMultipleStates_hh
#define COOLFluiD_Muffin_ComMultipleStates_hh

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/Node.hh"
#include "Framework/State.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Command for queuing multiple states
class ComMultipleStates : public MuffinCom {

 public:  // core functions

  /// Multiple states queue command constructor
  ComMultipleStates(const std::string& name);

  /// Multiple states queue command destructor
  ~ComMultipleStates() {}

  /// Defines the Config Option's of this command
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase (virtual implementation)
  void setup();

  /// Iterate on solution and derivatives queues (virtual implementation)
  void execute();


 private:  // intrinsic functions

  /// Calculate State derivatives, for a given TRS
  void setStatesDerivatives(Framework::TopologicalRegionSet& trs, Framework::DataHandle< CFreal >& ddx, Framework::DataHandle< CFreal >& ddy, Framework::DataHandle< CFreal >& ddz);


 private:  // sockets

  /// Socket to access nodes
  Framework::DataSocketSink < Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink < Framework::State*,Framework::GLOBAL > s_states;

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_nvolume;

  /// Sockets to provide solution at t-0 and (nodal) spacial derivative x, y and z components
  Framework::DataSocketSource< CFreal > s_tm0_sol;
  Framework::DataSocketSource< CFreal > s_tm0_ddx;
  Framework::DataSocketSource< CFreal > s_tm0_ddy;
  Framework::DataSocketSource< CFreal > s_tm0_ddz;

  /// Sockets to provide solution at t-1 and (nodal) spacial derivative x, y and z components
  Framework::DataSocketSource< CFreal > s_tm1_sol;
  Framework::DataSocketSource< CFreal > s_tm1_ddx;
  Framework::DataSocketSource< CFreal > s_tm1_ddy;
  Framework::DataSocketSource< CFreal > s_tm1_ddz;

  /// Sockets to provide solution at t-2 and (nodal) spacial derivative x, y and z components
  Framework::DataSocketSource< CFreal > s_tm2_sol;
  Framework::DataSocketSource< CFreal > s_tm2_ddx;
  Framework::DataSocketSource< CFreal > s_tm2_ddy;
  Framework::DataSocketSource< CFreal > s_tm2_ddz;


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
    r.push_back(&s_tm0_sol);
    r.push_back(&s_tm0_ddx);
    r.push_back(&s_tm0_ddy);
    r.push_back(&s_tm0_ddz);
    r.push_back(&s_tm1_sol);
    r.push_back(&s_tm1_ddx);
    r.push_back(&s_tm1_ddy);
    r.push_back(&s_tm1_ddz);
    r.push_back(&s_tm2_sol);
    r.push_back(&s_tm2_ddx);
    r.push_back(&s_tm2_ddy);
    r.push_back(&s_tm2_ddz);
    return r;
  }


 private:  // data

  /// Number of solutions to queue (configurable)
  CFuint m_qsolution;

  /// Number of solution spacial derivatives to queue (configurable)
  CFuint m_qderivatives;

  /// Number of nodes
  CFuint nbNodes;

  /// Number of states
  CFuint nbStates;

  /// Number of dimensions (Node size)
  CFuint nbDim;

  /// Number of variables per state (State size)
  CFuint nbVar;

  /// If this is the first ::execute call
  bool m_firstexecute;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_ComMultipleStates_hh
