#ifndef COOLFluiD_Numerics_FiniteVolume_StdALESetup_hh
#define COOLFluiD_Numerics_FiniteVolume_StdALESetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a Finite Volume Method using ALE formulation.
 * It complements the non-ALE Setup.
 *
 * @author Thomas Wuilbaut
 */
class StdALESetup : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit StdALESetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdALESetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Socket for the past coordinates of the Node's
  Framework::DataSocketSource<Framework::Node*> socket_pastNodes;

  /// Socket for the future coordinates of the Node's
  Framework::DataSocketSource<Framework::Node*> socket_futureNodes;

  /// Socket for the past volumes of cells
  Framework::DataSocketSource<CFreal> socket_pastVolumes;

  /// Socket for the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// Socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// Socket for the volumes of cells
  Framework::DataSocketSink<CFreal> socket_volumes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdALESetup_hh

