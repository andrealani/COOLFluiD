#ifndef COOLFluiD_UFEM_Solve_hh
#define COOLFluiD_UFEM_Solve_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "UFEM/UFEMSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a UFEM Command to discretize the inner cells.
class UFEM_API Solve : public UFEMSolverCom {
public:

  /// Constructor.
  explicit Solve(const std::string& name);

  /// Destructor.
  ~Solve();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Execute Processing actions
  void executeOnTrs();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// socket for rhs
  Framework::DataSocketSink< CFreal> socket_rhs;
  /// socket for states
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL > socket_states;
  Framework::DataSocketSink< Framework::State* > socket_interStates;

  /// number of equations
  CFuint m_nbEqs;
  /// number of dimensions
  CFuint m_dim;

}; // class Solve

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_Solve_hh

