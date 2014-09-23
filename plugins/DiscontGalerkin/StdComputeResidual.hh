#ifndef COOLFluiD_DiscontGalerkin_StdComputeResidual_hh
#define COOLFluiD_DiscontGalerkin_StdComputeResidual_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"

#include "DiscontGalerkin/DGElemTypeData.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This command to fill RHS after solution of linear system
 * @author Martin Holik
 */
class StdComputeResidual : public DiscontGalerkinSolverCom {
public: // functions

  /// Constructor
  explicit StdComputeResidual(const std::string& name);

  /// Destructor
  virtual ~StdComputeResidual();

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /// the socket for the norm computation
  Framework::DataSocketSource <CFreal>  socket_rhs_norm;
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State*, Framework::GLOBAL >  socket_states;
  /// the socket to the data handle of the old state's
  Framework::DataSocketSink < Framework::State* >  socket_old_states;

private: // data

}; // class StdComputeResidual

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_StdComputeResidual_hh

