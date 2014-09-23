#ifndef COOLFluiD_Numerics_FluctSplit_StdALEUpdate_hh
#define COOLFluiD_Numerics_FluctSplit_StdALEUpdate_hh

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
class StdALEUpdate : public FluctuationSplitCom {
public:

  /**
   * Constructor.
   */
  explicit StdALEUpdate(const std::string& name);

  /**
   * Destructor.
   */
  ~StdALEUpdate();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // helper method

  /**
   * Compute the speed of the cell
   */
  void computeCellSpeed();

  /**
   * Update the Normals
   */
  void updateNormalsData();

  /**
   * Store the Normals at past time
   */
  void savePastNormals();

protected:

  /// Socket to store the past normals
  Framework::DataSocketSink<
                              InwardNormalsData*> socket_normals;

  /// Socket to store the past normals
  Framework::DataSocketSink<
                              InwardNormalsData*> socket_pastNormals;

  /// Socket for info on the cell speed
  Framework::DataSocketSink<RealVector> socket_cellSpeed;

}; // class Update

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdALEUpdate_hh

