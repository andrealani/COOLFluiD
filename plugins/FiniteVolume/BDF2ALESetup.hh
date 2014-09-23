#ifndef COOLFluiD_Numerics_FiniteVolume_BDF2ALESetup_hh
#define COOLFluiD_Numerics_FiniteVolume_BDF2ALESetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdALESetup.hh"
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
class BDF2ALESetup : public StdALESetup {
public:

  /**
   * Constructor.
   */
  explicit BDF2ALESetup(const std::string& name);

  /**
   * Destructor.
   */
  ~BDF2ALESetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

protected:

  /// Socket for the past coordinates of the Node's (at time n-1)
  Framework::DataSocketSource<Framework::Node*> socket_pastPastNodes;

  /// Socket for the nodes
  Framework::DataSocketSource<CFreal> socket_avNormals;

  /// Socket for the Normals
  Framework::DataSocketSink<CFreal> socket_normals;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BDF2ALESetup_hh

