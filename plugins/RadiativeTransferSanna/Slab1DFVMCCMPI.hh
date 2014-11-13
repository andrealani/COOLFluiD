#ifndef COOLFluiD_Transition_Slab1DFVMCCMPI_hh
#define COOLFluiD_Transition_Slab1DFVMCCMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMultiMap.hh"

#include "RadiativeTransferSanna/Slab1DFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class RadiationLibrary;
  }
  
  namespace Numerics {
    namespace FiniteVolume {
      class CellCenterFVMData;
    }
  }
  
  namespace RadiativeTransfer {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 *
 * @author Andrea Lani
 * @author Alessandro Munafo'
 *
 */
class Slab1DFVMCCMPI : public Slab1DFVMCC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Slab1DFVMCCMPI(const std::string& name);

  /**
   * Default destructor
   */
  ~Slab1DFVMCCMPI();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();
  
protected:

  /**
   * Execute on a set of dofs
   */
  virtual void executeOnTrs();
   
  /**
   * Build mesh lines
   */
  virtual void buildMeshLines();
  
  /**
   * This class groups data defining a partition face
   *
   * @author Andrea Lani
   *
   */
  struct PartitionFaceData {
    std::vector<CFuint> nodeIDs;
    std::vector<CFuint> nbFaceNodes;
    std::vector<CFuint> globalCellIDs;
    std::vector<int> globalFaceIDs; // new entry
    std::vector<CFuint> nbLayersToAdd;
    
    const PartitionFaceData& operator= (const PartitionFaceData& d)
    {
      nodeIDs = d.nodeIDs;
      nbFaceNodes = d.nbFaceNodes;
      globalCellIDs = d.globalCellIDs;
      globalFaceIDs = d.globalFaceIDs; // new entry
      nbLayersToAdd = d.nbLayersToAdd;
      return *this;
    }
  };
  
  /**
    * Process the boundary faces
    */
  void processBFaces(const std::vector<std::string>& trsNames,
		     Common::CFMultiMap<CFuint,CFuint>& mapGlobalIDToLocalFaceID,
		     const std::vector<bool>& isPartitionFace,
		     std::vector<bool>& faceToConsider,
		     std::vector<CFuint>& nbLayersToAdd,
		     std::vector<CFint>& globalCellIDs,
		     std::vector<CFint>& globalFaceIDs,
		     std::vector<std::vector<CFint> >& lineGlobalCellIDs,
		     std::vector<RealVector>& distanceFromStagPoint);
  
  /**
   * Scan the B layers and build the corresponding face data
   */
  void scanLayers(Common::SafePtr<Framework::TopologicalRegionSet> wallTRS,
		  const std::vector<bool>& isPartitionFace,
		  const std::vector<bool>& faceToConsider,
		  const std::vector<CFuint>& nbLayersToAdd,
		  const std::vector<CFint>& globalCellIDs,
		  const std::vector<CFint>& globalFaceIDs,
		  std::vector<CFint>& lineGlobalCellIDsInTRS,
		  RealVector& distanceFromStagPoint,
		  PartitionFaceData& sendData);
  
private:
  
  typedef Common::CFMultiMap<CFuint,CFuint>::MapIterator MapItr; 
  
  /// rank of current processor
  CFuint m_myRank;
  
  /// local number of lines
  std::vector<CFuint> m_localNbLines;
  
  /// number of states in each processor
  std::vector<CFuint> m_nbStatesInProc;
    
  /// global state IDs in this processor
  std::vector<CFuint> m_sendGlobalIDs;
  
  /// temporary array for face midpoint
  RealVector m_midFace;
  
  /// map the TRS name to the number of faces on such a TRS
  Common::CFMap<std::string, CFuint> m_mapTrsNameToNbFaces;
  
  /// map the global cell (state) ID to a local mesh line ID
  Common::CFMultiMap<CFuint, CFuint> m_mapGlobalToMeshLineID;
     
  /// map global state IDs to local stateIDs in this process
  Common::CFMultiMap<CFuint, CFuint> m_mapGlobalToLocalStateIDs;
  
  /// local (inside the InnerFaces TRS) inner faceID to global wall face ID
  std::vector<int> m_innerFaceIDToGlobalWallFaceID;
    
  /// array with states and corrresponding nodes data
  std::vector<CFreal> m_statesNodes;
  
  /// send array
  std::vector<CFreal> m_sendArray;
  
  /// local qrad array corresponding to the llocal portion of the 
  /// structured mesh by line
  std::vector<CFreal> m_qRadByLine;
  
  /// qrad array to send
  std::vector<CFreal> m_qRadToSend;
  
  /// qrad array to receive
  std::vector<CFreal> m_qRadToRecv;
  
  /// array storing the destination socket_qrad (state) local IDs per rank
  std::vector<std::vector<CFuint> > m_destIDsPerRank;
  
  /// array storing the local IDs in the m_qRadByLine to be sent to each rank
  std::vector<std::vector<CFuint> > m_donorIDsPerRank;
  
  /// destination local IDs to be sent
  std::vector<CFint> m_destIDsToSend;
  
  /// destination local IDs to be received
  std::vector<CFint> m_destIDsToRecv;
  
  /// array storing the data count to send
  std::vector<int> m_sendCount;
  
  /// array storing the displacement of the data to send
  std::vector<int> m_sendDispl;
  
  /// array storing the data count to receive
  std::vector<int> m_recvCount;
  
  /// array storing the displacement of the data to receive
  std::vector<int> m_recvDispl;
  
  /// maximum number of normal faces in the direction orthogonal to the wall
  CFuint m_maxNbNormalFaces;
  
  /// maximum number of cells in the direction orthogonal to the wall
  CFuint m_nbCellsOnLine;
  
}; /// end of class Slab1DFVMCCMPI

//////////////////////////////////////////////////////////////////////////////

    } /// namespace RadiativeTransfer

} /// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /// COOLFluiD_RadiativeTransfer_Slab1DFVMCCMPI_hh
