#ifndef COOLFluiD_Framework_BaseSetupFVMCC_hh
#define COOLFluiD_Framework_BaseSetupFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ProxyDofIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a basic setup command creating data needed by a basic 
  * cell-centered Finite Volume solver
  *
  * @author Andrea Lani
  */
 template <typename BASE>
 class BaseSetupFVMCC : public BASE {
 public:

   /**
    * Defines the Config Option's of this class
    * @param options a OptionList where to add the Option's
    */
   static void defineConfigOptions(Config::OptionList& options);

   /**
    * Constructor.
    */
   explicit BaseSetupFVMCC(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BaseSetupFVMCC();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );
   
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
   
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
   * Allocate arrays
   */
  virtual void allocateArrays();
   
  /**
   * Compute the normals data
   */
  virtual void computeNormalsData();

  /**
   * Compute the cell volumes
   */
  virtual void computeCellVolumes();

  /**
   * Store the (outward) normals at the boundaries
   */
  virtual void storeBoundaryNormals();
   
  /// assign global ghost state IDs (corresponding to existing cell-state IDs) for partition faces 
  virtual void assignPartitionFaceGlobalGhostStateIDS();
   
protected:
  
  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_dynamicSockets;
   
  /// storage of face normals
  Framework::DataSocketSource<CFreal> socket_normals;
  
  /// storage of face areas
  Framework::DataSocketSource<CFreal> socket_faceAreas;
  
  /// storage of face centroids
  Framework::DataSocketSource<CFreal> socket_faceCenters;
  
  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSource<CFint> socket_isOutward;
  
  /// flags for the cells
  Framework::DataSocketSource<bool> socket_cellFlag;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSource<Framework::State*> socket_gstates;
  
  /// storage of the cell volumes
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
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for delayed configuration of an object (MeshFormatConverter)
  Config::ConfigArgs m_stored_args;
  
}; // class BaseSetupFVMCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseSetupFVMCC.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseSetupFVMCC_hh

