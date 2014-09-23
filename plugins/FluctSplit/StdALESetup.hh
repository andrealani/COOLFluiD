#ifndef COOLFluiD_Numerics_FluctSplit_StdALESetup_hh
#define COOLFluiD_Numerics_FluctSplit_StdALESetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to setup the MeshData.
 *
 * @author Thomas Wuilbaut
 */
class StdALESetup : public FluctuationSplitCom {
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
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

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
  Framework::DataSocketSource<
                              InwardNormalsData*> socket_pastNormals;

  /// Socket to store the past normals data
  Framework::DataSocketSource<
                              CFreal> socket_pastNormalsData;

  /// Socket for info on the cell past volume
  Framework::DataSocketSource<CFreal> socket_pastCellVolume;

  /// Socket for info on the cell speed
  Framework::DataSocketSource<RealVector> socket_cellSpeed;

  /// Socket with past coordinates of the nodes
  /// which are a copy of the original nodes
  Framework::DataSocketSource<
                              Framework::Node*> socket_pastNodes;

  // the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  // the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdALESetup_hh

