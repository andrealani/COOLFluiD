#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferFVSolar_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferFVSolar_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Environment {
    class FileHandlerInput;
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace RadiativeTransfer {
    class RadiationPhysicsHandler;
    
//////////////////////////////////////////////////////////////////////////

/**
 * This class compute the radiative heat transfer using a Finite Volume algorithm
 *
 * @author Andrea Lani 
 * @author Michaela Brchnelova 
 */
class RadiativeTransferFVSolar : public Framework::DataProcessingCom {
public:
  
  /**
   * Constructor
   */
  RadiativeTransferFVSolar(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RadiativeTransferFVSolar();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );
  
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
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();
  
  /**
   * Execute on a set of dofs
   */
  virtual void execute();
  
  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected: //function

  /// store the cell IDs of the paths for integration 
  void storeIntegralPathIDs();
    
  /// set the normal corresponding to the given face ID
  void setFaceNormal(const CFuint faceID, const CFuint elemID) 
  {
    Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle();
    const CFuint startID = faceID*3;
    const CFreal factor = (static_cast<CFuint>
			   (socket_isOutward.getDataHandle()[faceID]) != elemID) ? -1. : 1.;
    for (CFuint dir = 0; dir < 3; ++dir) {
      m_normal[dir] = normals[startID+dir]*factor;
    }    
  }
    
  /// get the neighbor cell ID to the given face and cell
  CFuint getNeighborCellID(const CFuint faceID, const CFuint cellID) const 
  {
    // find the TRS to which the current face belong
    const Framework::TopologicalRegionSet& faceTrs = *m_mapGeoToTrs->getTrs(faceID);
    // find the local index for such a face within the corresponding TRS
    const CFuint faceIdx = m_mapGeoToTrs->getIdxInTrs(faceID);
    // first neighbor of the current face
    const CFuint sID0 = faceTrs.getStateID(faceIdx, 0);
    return (cellID == sID0) ? faceTrs.getStateID(faceIdx, 1) : sID0;
  }
    
  /// @return the wall face ID to be used inside the qradFluxWall socket
  /// @post returns -1 if the face is not a wall face
  CFint getWallFaceID(const CFuint faceID)
  {
    if (m_isWallFace[faceID]) {
      // find the local index for such a face within the corresponding TRS
      const CFuint faceIdx = m_mapGeoToTrs->getIdxInTrs(faceID);
      const std::string wallTRSName = m_mapGeoToTrs->getTrs(faceID)->getName();
      return m_mapWallTRSToOffsetFaceID.find(wallTRSName)->second + faceIdx;
    }
    return -1;
  }

  /// Gets the opposite local face ID (in 2D/3D structured meshes or boundary layers)
  CFuint getOppositeIFace(CFuint iFace, CFuint dim, CFuint nbCellNodes) const;
  
protected: //data
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;  
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;  
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;  
  
  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward; 
  
  /// storage of normals
  Framework::DataSocketSink<CFreal> socket_normals; 
  
  /// storage of face centers
  Framework::DataSocketSink<CFreal> socket_faceCenters; 
  
  /// storage of face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas; 
  
  /// storage of the divq 
  Framework::DataSocketSource <CFreal> socket_divq;
  
  /// storage of the qx 
  Framework::DataSocketSource <CFreal> socket_qx;
  
  /// storage of the qy
  Framework::DataSocketSource <CFreal> socket_qy;
  
  /// storage of the qz
  Framework::DataSocketSource <CFreal> socket_qz;
  
  /// the socket to the radiative heat flux at the wall faces
  Framework::DataSocketSource <CFreal> socket_qradFluxWall;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library; 
  
  /// pointer to the radiation library interface
  Common::SharedPtr<RadiationPhysicsHandler> m_radiation;
  
  /// map faces to corresponding TRS and index inside that TRS
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs;
  
  /// flag array telling whether a face is on the wall
  std::vector<bool> m_isWallFace;

  /// array storing the wall face ID corresponding to the given cellID
  std::vector<CFuint> m_cellID2WallFaceID;
  
  /// table storing cellIDs belonging to grid lines starting in each wall faceID
  Common::Table<CFuint> m_tableWallFaceID2CellIDs;
  
  /// map with TRS name as key and the offset for the wall face IDs as value
  std::map<std::string, CFuint> m_mapWallTRSToOffsetFaceID;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// wall face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_wallFaceBuilder;

  /// temporary normal to the face
  RealVector m_normal; 

  /// maximum number of normal faces to neglect for carbuncle fix
  CFuint m_maxNbNormalFaces;
  
  /// names of the TRSs of type "Wall"
  std::vector<std::string> m_wallTrsNames;
  
  /// name of the radiation namespace
  std::string m_radNamespace;
  
}; // end of class RadiativeTransferFVSolar
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_RadiativeTransferFVSolar_hh



