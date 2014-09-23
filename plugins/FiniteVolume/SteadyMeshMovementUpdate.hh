#ifndef COOLFluiD_Numerics_FiniteVolume_SteadyMeshMovementUpdate_hh
#define COOLFluiD_Numerics_FiniteVolume_SteadyMeshMovementUpdate_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command to be executed after
   * the mesh has been updated
   */

//////////////////////////////////////////////////////////////////////////////

class SteadyMeshMovementUpdate : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit SteadyMeshMovementUpdate(const std::string& name);

  /**
   * Destructor.
   */
  ~SteadyMeshMovementUpdate()
  {
  }

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

  /**
   * Compute the nodes at intermediate time
   * for the cell center and ghost states
   */
  void modifyOffMeshNodes();

  /**
   * Backup the face area at time n
   */
  void backupFaceAreas();

  /**
   * Reset IsOutward to -1
   */
  void resetIsOutward();

  /**
   * Update the normals (at intermediate time)
   */
  virtual void updateNormalsData();

  /**
   * Update the faceAreas (at intermediate time)
   */
  void updateFaceAreas();

  /**
   * Compute the cell volumes
   */
  void updateCellVolume();

  /**
   * Update the reconstructor...
   */
  void updateReconstructor();

protected:

  // handle to nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // handle to the flag telling the ID of the cell for which the normal is outward
  Framework::DataSocketSink< CFint> socket_isOutward;

  // handle to the normals
  Framework::DataSocketSink< CFreal> socket_normals;
 
  // handle to the face areas
  Framework::DataSocketSink< CFreal> socket_faceAreas;
  
  // handle to the volumes
  Framework::DataSocketSink< CFreal> socket_volumes;

  // handle to the ghost states
  Framework::DataSocketSink< Framework::State*> socket_gstates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SteadyMeshMovementUpdate_hh

