#ifndef COOLFluiD_Numerics_FluctSplit_StdALEUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_StdALEUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Thomas Wuilbaut
 */
class StdALEUnSetup : public FluctuationSplitCom {
public:

  /**
   * Constructor.
   */
  explicit StdALEUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdALEUnSetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected:

  /// Socket to store the past normals
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_pastNormals;

  /// Socket with past coordinates of the nodes
  /// which are a copy of the original nodes
  Framework::DataSocketSink<
                            Framework::Node*> socket_pastNodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdALEUnSetup_hh

