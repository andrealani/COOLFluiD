#ifndef COOLFluiD_FluxReconstructionCUDA_StencilCUDASetup_hh
#define COOLFluiD_FluxReconstructionCUDA_StencilCUDASetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ProxyDofIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a standard FR Method when CUDA is used.
 *
 * @author Andrea Lani
 */
class StencilCUDASetup : public FluxReconstructionSolverCom {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit StencilCUDASetup(const std::string& name);
  
  /**
   * Destructor.
   */
  ~StencilCUDASetup();
  
  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // helper method

  /**
   * Compute the normals data
   */
  void computeNormalsData();
  
  /// assign global ghost state IDs (corresponding to existing cell-state IDs) for partition faces 
  void assignPartitionFaceGlobalGhostStateIDS();
  
  /**
   * Compute and store the stencil for each cell
   */
  void computeStencil();
  
protected:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _dynamicSockets;
  
  /// flags for the cells (to be removed)
  Framework::DataSocketSource<bool> socket_cellFlag;
  
  /// storage of face areas (to be removed)
  Framework::DataSocketSource<CFreal> socket_faceAreas;
  
  /// @todo missing documentation
  Framework::DataSocketSource<CFint> socket_trsID;

  /// rank of the processor containing the right state for partition faces
  Framework::DataSocketSource<CFuint> socket_rankPartitionFaces;
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSource<std::vector<Framework::State*> > socket_stencil;
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for delayed configuration of an object (MeshFormatConverter)
  Config::ConfigArgs m_stored_args;
  
  // type of the stencil
  std::string _stencilType;
      
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionCUDA_StencilCUDASetup_hh
