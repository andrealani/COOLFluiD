#ifndef COOLFluiD_Numerics_FiniteVolume_StdLinSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_StdLinSetup_hh

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
class StdLinSetup : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit StdLinSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdLinSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected:

  /// Socket for the ghost states
  Framework::DataSocketSource<Framework::State*> socket_linearizedGhostStates;

  /// Socket for the past ghost states
  Framework::DataSocketSource<Framework::State*> socket_pastGhostStates;

  /// Socket for the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdLinSetup_hh

