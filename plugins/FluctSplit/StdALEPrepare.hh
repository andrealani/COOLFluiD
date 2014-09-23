#ifndef COOLFluiD_Numerics_FluctSplit_StdALEPrepare_hh
#define COOLFluiD_Numerics_FluctSplit_StdALEPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to Update the MeshData.
   */
class StdALEPrepare : public FluctuationSplitCom {
public:

  /**
   * Constructor.
   */
  explicit StdALEPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~StdALEPrepare()
  {
  }

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private: // helper method

  /**
   * Update the Normals
   */
  void updateNormalsData();

  /**
   * Store the Normals at past time
   */
  void savePastNormals();

private:

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;
  
  /// socket for pastNode's
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;
  
  /// socket for pastNode's
  Framework::DataSocketSink<CFreal> socket_pastCellVolume;
  
  /// socket for pastNode's
  Framework::DataSocketSink<RealVector> socket_cellSpeed;

}; // class Update

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdALEPrepare_hh

