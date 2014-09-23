#ifndef COOLFluiD_Muffin_BC_hh
#define COOLFluiD_Muffin_BC_hh

#include "Common/COOLFluiD.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Base class for boundary condition
class BC : public MuffinCom,
           public Common::TaggedObject {

 public:  // functions

  /// Boundary condition constructor
  BC(const std::string& name);

  /// Boundary condition destructor
  ~BC();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  virtual void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  virtual void setup();

  /// Apply boundary condition (wrapper for appropriate applyOnSystem...)
  virtual void apply(const Common::SafePtr< System > s);

  /// Apply boundary condition (specific system application, to be overridden)
  virtual void applyOnSystemFlow(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t) {}
  virtual void applyOnSystemTemp(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t) {}
  virtual void applyOnSystemTurb(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t) {}
  virtual void applyOnSystemMITReM(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t) {}


 protected:  // functions

  /// Log information, debug and error messages
  void log(const std::string& msg) { getMethodData().log("BC " + getName() + ": " + msg); }
  void ver(const std::string& msg) { getMethodData().ver("BC " + getName() + ": " + msg); }
  void err(const std::string& msg) { getMethodData().err("BC " + getName() + ": " + msg); }

  /// Set var_type states at node, to initialize the solution field
  void setInitialState(const CFuint& n, const var_type& state, const CFreal& v);

  /// Set var_type states at node, to initialize the solution field
  void setInitialState(const CFuint& n, const var_type& state, const RealVector& vv);

  /// Set linear system equation as Dirichlet condition (with equation factor)
  void setDirichletCondition( System& s,
    const CFuint n, const CFuint e, const double v,
    const double f=1. );

  /// Set linear system equation as Dirichlet condition, coupling node 1,
  /// equation 1 to node 2, equation 2 (with equation factor)
  void setDirichletCondition( System& s,
    const CFuint n1, const CFuint e1, const double v1,
    const CFuint n2, const CFuint e2, const double v2,
    const double f=1. );


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
    r.push_back(&s_mn_priority);
    return r;
  }

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets()
  {
    return std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >();
  }


 protected:  // sockets

  /// Socket to access RHS
  Framework::DataSocketSink< CFreal > s_rhs;

  /// Socket to access nodes
  Framework::DataSocketSink< Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink< Framework::State*,Framework::GLOBAL > s_states;

  /// Socket to access mapping surface face to corresponding boundary element
  Framework::DataSocketSink< std::pair< CFuint,CFuint > > s_faceneighcell;

  /// Socket to access node-wise boundary condition application priority
  Framework::DataSocketSink< CFuint > s_mn_priority;


 public:  // data

  /// Priority in applicable nodes over other boundary conditions
  CFuint m_bnpriority;

  /// Apply to a single point coordinates
  std::vector< double > m_point_xyz;

  /// Apply to a single point node index
  int m_point_index;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BC_hh

