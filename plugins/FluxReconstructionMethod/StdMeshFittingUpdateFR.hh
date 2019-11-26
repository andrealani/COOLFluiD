#ifndef COOLFluiD_FluxReconstructionMethod_StdMeshFittingUpdateFR_hh
#define COOLFluiD_FluxReconstructionMethod_StdMeshFittingUpdateFR_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolver.hh" 
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command to be executed after
   * the mesh has been updated
   */

//////////////////////////////////////////////////////////////////////////////

class StdMeshFittingUpdateFR : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit StdMeshFittingUpdateFR(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdMeshFittingUpdateFR();
  
  /**
   * Execute Processing actions
   */
  virtual void execute();
  /**
   * Creates socket with the cell volume for each state
   */
  void computeStatesVolumes();


  /**
   * Computes the unit normals and the the size of the projection vectors in face flux points.
   */
  void computeFaceJacobianVectorSizes();
   /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */


  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data


  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for cell volumes
  Framework::DataSocketSink< CFreal > socket_volumes;

  /// socket for size of face normal jacobian in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;
  
  /// socket for normals
  Framework::DataSocketSink< CFreal > socket_normals;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdMeshFittingUpdate_hh
