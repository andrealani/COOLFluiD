#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingCorrFluxMPI_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingCorrFluxMPI_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/CFMultiMap.hh"
#include "FiniteVolumeNavierStokes/StegerWarmingCorrFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class Node;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the StegerWarmingCorrMPI flux
 *
 * @author Andrea Lani
 *
 */
class StegerWarmingCorrFluxMPI : public StegerWarmingCorrFlux {
public:
  
  struct PartitionFaceData {
    std::vector<CFuint> nodeIDs;
    std::vector<CFuint> nbFaceNodes;
    std::vector<CFuint> globalCellIDs;
    std::vector<CFuint> nbLayersToAdd;
    
    const PartitionFaceData& operator= (const PartitionFaceData& d)
    {
      nodeIDs = d.nodeIDs;
      nbFaceNodes = d.nbFaceNodes;
      globalCellIDs = d.globalCellIDs;
      nbLayersToAdd = d.nbLayersToAdd;
      return *this;
    }
  };
  
  /**
   * Constructor
   */
  StegerWarmingCorrFluxMPI(const std::string& name);

  /**
   * Default destructor
   */
  ~StegerWarmingCorrFluxMPI();

  /**
   * Build the face data
   */
  virtual void buildFaceBCData();

private:
  
  /**
    * Process the boundary faces
    */
  void processBFaces(const std::vector<std::string>& trsNames,
		     Common::CFMultiMap<CFuint,CFuint>& mapGlobalIDToLocalFaceID,
		     const std::vector<bool>& isPartitionFace,
		     std::vector<bool>& faceToConsider,
		     std::vector<CFuint>& nbLayersToAdd,
		     std::vector<CFint>& globalCellIDs);
  
  /**
   * Scan the B layers and build the corresponding face data
   */
  void scanLayers(Common::SafePtr<Framework::TopologicalRegionSet> wallTRS,
		  const std::vector<bool>& isPartitionFace,
		  const std::vector<bool>& faceToConsider,
		  const std::vector<CFuint>& nbLayersToAdd,
		  const std::vector<CFint>& globalCellIDs,
		  PartitionFaceData& sendData);
  
private:
    
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> _cellBuilder;
  
  
}; // end of class StegerWarmingCorrFluxMPI

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingCorrFluxMPI_hh
