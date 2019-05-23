#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRHSCell_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRHSCell_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRHS.hh"
#include "FiniteVolume/KernelData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class CellConn;
  }
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes with CUDA bindings
 *
 * @author Andrea Lani
 *
 */

//This template include one more parameter, the source term
template <typename SCHEME, typename PHYSICS, typename SOURCE, typename POLYREC, typename LIMITER, CFuint NB_BLOCK_THREADS>
class FVMCC_ComputeSourceRHSCell: public FVMCC_ComputeRHS {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeSourceRHSCell(const std::string& name);
  
  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeSourceRHSCell();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Un Setup private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
      
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:
  
  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Store the stencil data
  virtual void storeStencilData();
  
  /// Copy the local connectivity data to GPU
  void copyLocalCellConnectivity();
  
protected:
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// socket for uZ values
  Framework::DataSocketSink<CFreal> socket_uZ;
 
  
  Framework::DataSocketSink<CFreal> socket_volumes;
  


  /// cell-face connectivity
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellFaces;
  
  /// cell-nodes connectivity
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellNodes;
  
  /// storage of the cell centers (AL: temporary solution) 
  CudaEnv::CudaVector<CFreal> m_centerNodes;
  
  /// storage of the ghost states (AL: temporary solution) 
  CudaEnv::CudaVector<CFreal> m_ghostStates;
  
  /// storage of the ghost nodes (AL: temporary solution) 
  CudaEnv::CudaVector<CFreal> m_ghostNodes;
  
  /// storage of useful cell info: 
  /// in [cellID*5+0] - ptr to corresponding stencil
  /// in [cellID*5+1] - stencil size
  /// in [cellID*5+2] - number of cell faces
  /// in [cellID*5+3] - cell geoshape
  /// in [cellID*5+4] - number of active cell faces (partition faces are excluded)  
  CudaEnv::CudaVector<CFuint, CudaEnv::MallocHostAlloc<CFuint> > m_cellInfo;
  
  /// stencil connectivity for cellID: 
  /// starts at m_cellInfo[cellID*5]
  /// its size is given by m_cellInfo[cellID*5+1]
  /// first m_cellInfo[cellID*5+2] are faces
  CudaEnv::CudaVector<CFuint, CudaEnv::MallocHostAlloc<CFuint> > m_cellStencil;
  
  /// storage of flags for neighbors (1: internal, 0:partition, <0: boundary)
  /// starts at m_cellInfo[cellID*5]
  /// its size is given by m_cellInfo[cellID*5+1]
  /// first m_cellInfo[cellID*5+2] are faces
  CudaEnv::CudaVector<CFint, CudaEnv::MallocHostAlloc<CFint> > m_neighborTypes;
  
  // cell connectivity
  CudaEnv::CudaVector<Framework::CellConn, CudaEnv::MallocHostAlloc<Framework::CellConn> > m_cellConn;
  
  /// number of blocks in x 
  CFuint m_nbBlocksPerGridX;
  
  /// number of blocks in y
  CFuint m_nbBlocksPerGridY;
  
  /// number of cells per block
  CFuint m_nbCellsPerBlock;
  
  /// flag telling to solve on GPU
  bool m_onGPU;

  /// build the mesh on the GPU for Paralution
  bool m_useParalutionPtr;  

}; // class FVMCC_ComputeSourceRHSCell

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeSourceRHSCell.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRHSCell_hh
