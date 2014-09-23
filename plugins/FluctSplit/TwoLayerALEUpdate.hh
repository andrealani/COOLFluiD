#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerALEUpdate_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerALEUpdate_hh

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
   *
   * @author Thomas Wuilbaut
   *
   */

class TwoLayerALEUpdate : public FluctuationSplitCom {
public:

  /**
   * Constructor.
   */
  explicit TwoLayerALEUpdate(const std::string& name);

  /**
   * Destructor.
   */
  ~TwoLayerALEUpdate();

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
   * Update the Normals
   */
  void updateNormalsData();

  /**
   * Store the Normals at past time
   */
  void savePastNormals();

  /**
   * Compute the Normals at intermediate time
   */
  void computeInterNormalsData();

  /**
   * Compute the Normals at intermediate time
   */
  void computeCellSpeed();

protected:

  /// Socket to store the past normals
  Framework::DataSocketSink<
                              InwardNormalsData*> socket_normals;

  /// Socket to store the past normals
  Framework::DataSocketSink<
                              InwardNormalsData*> socket_pastNormals;

  /// Socket to store the past normals
  Framework::DataSocketSink<
                              InwardNormalsData*> socket_interNormals;

  /// Socket for info on the cell speed
  Framework::DataSocketSink<RealVector> socket_cellSpeed;


}; // class Update

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerALEUpdate_hh

