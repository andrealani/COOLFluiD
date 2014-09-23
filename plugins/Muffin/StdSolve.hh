#ifndef COOLFluiD_Muffin_StdSolve_hh
#define COOLFluiD_Muffin_StdSolve_hh

#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// This is a standard command to solve a testcase using Muffin
/// (it simply executes the first Loop)
class StdSolve : public MuffinCom {

public:  // members

  /// Constructor
  explicit StdSolve(const std::string& name) :
      MuffinCom(name),
      s_rhs("rhs"),
      s_states("states")
  {}

  /// Destructor
  ~StdSolve()
  {}

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_rhs);
    r.push_back(&s_states);
    return r;
  }


 private:  // data

  /// Socket to access  RHS
  Framework::DataSocketSink< CFreal > s_rhs;

  /// Socket to access  states
  Framework::DataSocketSink< Framework::State*,Framework::GLOBAL > s_states;

}; // class Solve


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_StdSolve_hh

