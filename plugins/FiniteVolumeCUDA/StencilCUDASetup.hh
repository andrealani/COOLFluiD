#ifndef COOLFluiD_Numerics_FiniteVolume_StencilCUDASetup_hh
#define COOLFluiD_Numerics_FiniteVolume_StencilCUDASetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ProxyDofIterator.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a standard Finite Volume Method when CUDA is used.
 *
 * @author Andrea Lani
 */
class StencilCUDASetup : public CellCenterFVMCom {
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
  
  /**
   * Compute the cell volumes
   */
  void computeCellVolumes();
  
  /// assign global ghost state IDs (corresponding to existing cell-state IDs) for partition faces 
  void assignPartitionFaceGlobalGhostStateIDS();
  
  /**
   * Compute and store the stencil for each cell
   */
  void computeStencil();
  
protected:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _dynamicSockets;
  
  /// IDs corresponding to the cell for which the normal point outward (to be removed)
  Framework::DataSocketSource<CFint> socket_isOutward;
  
  /// flags for the cells (to be removed)
  Framework::DataSocketSource<bool> socket_cellFlag;
  
  /// storage of face normals (to be removed)
  Framework::DataSocketSource<CFreal> socket_normals;
  
  /// storage of face areas (to be removed)
  Framework::DataSocketSource<CFreal> socket_faceAreas;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSource<Framework::State*> socket_gstates;
  
  /// storage of the cell volumes (to be removed)
  Framework::DataSocketSource<CFreal> socket_volumes;
  
  /// Storage of the limiters
  /// there is one for each cell
  Framework::DataSocketSource<CFreal> socket_limiter;
  
  /// storage for interpolated values in the mesh vertices
  Framework::DataSocketSource<RealVector> socket_nstates;
  
  /// socket for proxy to be able to use the data handle of nodal states uniformly
  /// independently from the actual storage type being RealVector or State*
  Framework::DataSocketSource<Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;
  
  /// @todo missing documentation
  Framework::DataSocketSource<CFint> socket_trsID;

  /// rank of the processor containing the right state for partition faces
  Framework::DataSocketSource<CFuint> socket_rankPartitionFaces;
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSource<std::vector<Framework::State*> > socket_stencil;

  /// flags telling if diffusion terms have to be computed for a given cell
  Framework::DataSocketSource<CFreal> socket_activeDiffusion;
  
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
  
  /// initial limiter socket name
  std::string _limiterSocketName;
      
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StencilCUDASetup_hh
