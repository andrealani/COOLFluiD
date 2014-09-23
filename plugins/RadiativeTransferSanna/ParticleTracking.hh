#ifndef COOLFluiD_RadiativeTransfer_ParticleTracking_hh
#define COOLFluiD_RadiativeTransfer_ParticleTracking_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>

#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

  struct GeoEntOut
  {
    CFint exitFaceID;
    CFint exitCellID;
  };
  /**
   * This struct models a ray
   *
   * @author Alessandro Sanna
   * @author Andrea Lani
   */
   struct Ray
   {
     CFreal direction[3];
     CFreal startPoint[3];
     CFreal KS;
     CFreal actualKS;
     CFuint startEntityID; // first cell attached to partition face
     CFuint emittingEntityID;
     CFint emittingEntityType;  // before: 1 if it is a cell, 0 if it is a face
                                // now: -1 for a cell; globalID > 0 for a face
     CFuint emittingProcessorRank;
     CFreal tt;
     CFreal wavelength;
   };


  /**
   * This class implements the algorithms for the computation of non local parameters
   *
   * @author Alessandro Sanna
   *
   */
class ParticleTracking
{
public:
  
  /// enumerators that define the information type
  enum InfoType {ENDPOINT =0, DIRECTION =1};
  
  /**
   * Constructor
   */
  ParticleTracking();
  
  /**
   * Default destructor
   */
   ~ParticleTracking();

  /**
   * build cell IDs map
   */
   void buildCellIDmap();


   GeoEntOut myAxiRayTracing(Ray &ray, CFint currentCellID);

  /**
   * setup start point (if it is inside a cell): use the barycentre of the start cell
   */
  void setupStartPoint(CFuint startCellID);
  
  /**
   * setup start point (if it is inside a cell): use some certain coordinates given in startPoint
   */
  void setupStartPoint(CFuint startCellID, CFreal* startPoint)
  {
    _startCellID = startCellID;
    // use the coordinate stored in the vector startPoint to set the start Point coordinates
    for (CFuint d = 0; d < _dim; ++d) {
      _startPoint[d] = startPoint[d];
    }
  }
  
  /**
   * setup start point (if it is on a face): use the barycentre of the start face
   */
  void setupStartPoint(CFuint startFaceID, const std::string& faceTrs);
  
  /**
   * case1: pathInformation = endPointCoordinate; InformationType = ENDPOINT;
   *
   *         use this case1 if you want ALL and ONLY the cells that are crossed 
   *         by a straight line that connect the start point with the end point.
   *
   *
   * case2: pathInformation = direction cosines; InformationType = DIRECTION;
   *
   *         use this case2 if you want to know the cell that:
   *                                                          1) is a neighbour of the cell containing the start point;
   *                                                          2) is crossed by a straight line passing through the start point with the specified direction.
   */
  void setupPath(CFreal* pathInformation, InfoType InformationType);

  /**
   * set the DataSockets
   */
   void setDataSockets(Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
                       Framework::DataSocketSink<Framework::State*> gstatesSocket,
                       Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> nodesSocket,
                       Framework::DataSocketSink < RealVector > nstatesSocket,
                       Framework::DataSocketSink< CFreal> normalsSocket,
                       Framework::DataSocketSink<CFint> isOutwardsSocket
                       );
  
  /**
   * tracking for: set=1
   */
  void tracking(std::vector<CFuint>& CellList, CFuint actualCellID);

  /**
   * tracking for: set=2
   */
  CFuint tracking(CFuint actualCellID)
  {
    //if(_set == 2){
    if(_dim == 2) { 
      return tracking2D(actualCellID);
    }
    return tracking3D(actualCellID);
    //} 
  }
  
  /**
   * tracking 2D
   */
   CFuint tracking2D(CFuint cellID);

  /**
   * tracking 3D
   */
   CFuint tracking3D(CFuint cellID);

  /**
   * sort nodes in 2D
   */
  void sortNodes2D(std::vector<CFreal>& NodeXcoordinate, std::vector<CFreal>& NodeYcoordinate, const RealVector& CellBarycentre, CFuint nFaces);

  /**
   * sort nodes in 3D
   */
  void sortNodes3D(std::vector<CFreal>& nodeXcoordinate, std::vector<CFreal>& nodeYcoordinate, std::vector<CFreal>& nodeZcoordinate, 
		   const RealVector& barycentreActualCell, const RealVector& barycentreActualFace, CFuint p, CFuint r);
  
  /**
   * T2L test
   */
  void T2Ltest(const std::vector<CFreal>& nodeXcoordinate, 
	       const std::vector<CFreal>& nodeYcoordinate, 
	       CFuint nFaces, std::vector<CFuint>& intersectedFaces)
  {
    intersectedFaces.resize(0, 0);
    for(CFuint i=0; i<nFaces; ++i){
      const CFreal L_first = (nodeXcoordinate[2*i] - _startPoint[0])*(_endPoint[1] - _startPoint[1]) -
	(nodeYcoordinate[2*i] - _startPoint[1])*(_endPoint[0] - _startPoint[0]);
      const CFreal L_second = (nodeXcoordinate[2*i+1] - _startPoint[0])*(_endPoint[1] - _startPoint[1]) -
	(nodeYcoordinate[2*i+1] - _startPoint[1])*(_endPoint[0] - _startPoint[0]);
      if( (L_first/L_second) < 0 ){
	intersectedFaces.push_back(i);
      }
    }
  }
  
  /**
   * P2L test
   */
  bool P2Ltest(const std::vector<CFreal>& nodeXcoordinate, 
	       const std::vector<CFreal>& nodeYcoordinate, 
	       const std::vector<CFuint>& intersectedFaces, 
	       CFuint & exitFace);
  
  /**
   * T2I test
   */
  bool T2Itest(const std::vector<CFreal>& nodeXcoordinate, 
	       const std::vector<CFreal>& nodeYcoordinate, 
	       const std::vector<CFreal>& nodeZcoordinate);
  
  /**
   * P2I test
   */
  bool P2Itest(const std::vector<CFreal>& nodeXcoordinate, 
	       const std::vector<CFreal>& nodeYcoordinate, 
	       const std::vector<CFreal>& nodeZcoordinate);
  
  /**
   * get start Point coordinate
   */
  const RealVector& getStartPoint() {return _startPoint;}
  
  /**
   * get end Point coordinate
   */
  const RealVector& getEndPoint() {return _endPoint;}
  
  /**
   * get internal Point coordinate
   */
  const RealVector& getInternalPoint() {return _internalPoint;}
  
  /**
   * Get exit face ID (<0 is not found)
   */
  CFint getExitFaceID() {return _ExitFaceID;}
  
  /**
   * get exit face geometric entity
   */
  CFuint getExitFaceGlobalGhostID() 
  {
    using namespace COOLFluiD::Framework;
    
    _cellBuilder.getDataGE().idx = _cellIdx;
    GeometricEntity *const cell = _cellBuilder.buildGE();

    cf_assert(_FaceNb < cell->nbNeighborGeos());
    GeometricEntity* ExitFace = cell->getNeighborGeo(_FaceNb);
    cf_assert(ExitFace->getState(1)->isGhost());
    CFLog(DEBUG_MAX, "ParticleTracking::getExitFaceGlobalGhostID() => before globalID\n");
    const CFuint ghostGlobalID = ExitFace->getState(1)->getGlobalID();  
    CFLog(DEBUG_MAX, "ParticleTracking::getExitFaceGlobalGhostID() => after globalID\n\n");
    _cellBuilder.releaseGE();

    return ghostGlobalID;
  }

   /**
   * get startCellID
   */
  CFuint getStartCellID() {return _startCellID;}
  
  /**
   * get wallFaceID
   */
  CFuint getWallFaceID() {return _wallFaceID;}
  
  /**
   * get the intersection point between the exitFace and the path of the particle
   */
  RealVector& getIntersectionPoint();

  CFuint  getExitFaceGlobalGhostID(GeoEntOut out);

  
  /// average point
  template <typename RTYPE>
  void computeAverage(const std::vector<Framework::Node*>& nodes, 
		      CFuint nbNodes,  RTYPE& result)
  {
    const CFreal invNbNodes = 1./static_cast<CFreal>(nbNodes);
    result = 0.;
    for(CFuint n=0; n<nbNodes; ++n){
      result += *(nodes[n])*invNbNodes;
    }
  } 
  
  /// global to local cellID
  CFuint getLocalCellID(CFuint globalCellID) {return _CellIDmap.find(globalCellID);}


  
private: //data
  
  /// storage of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of the nodal state's
  Framework::DataSocketSink < RealVector > socket_nstates;

  /// handle to the face normals
  Framework::DataSocketSink< CFreal> socket_normals;

  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFint> socket_isOutward;
 
  /// cell builder
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> _cellBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
  
  /// std::vector defining the direction of the particle
  RealVector _direction;

  /// std::vector storing the cartesian coordinate of the ending point of the particle
  RealVector _endPoint;
  
  /// std::vector storing the cartesian coordinate of the starting point of the particle
  RealVector _startPoint;
  
  /// std::vector storing the cartesian coordinate of the internal point of the particle
  RealVector _internalPoint;
  
  /// std::vector storing the intersection point between the exitFace and the path of the particle
  RealVector _intersectionPoint;
  
  /// std::vector storing the baricentre of the cell
  RealVector _barycentreActualCell;
  
  /// std::vector storing the baricentre of the face
  RealVector _barycentreActualFace;
  
  /// arrays 
  RealVector _Lx;
  RealVector _Ly;
  RealVector _Lz;
  RealVector _Rx;
  RealVector _Ry;
  RealVector _Rz;
  RealVector _nn;
  RealVector _barycentreFace;
  
  std::vector<CFreal> _nodeXcoordinate;
  std::vector<CFreal> _nodeYcoordinate;
  std::vector<CFreal> _nodeZcoordinate;
  
  /// says true if the cell containing the exitPoint has been found (used if _set = 1)
  bool _foundIt;
  
  /// number of dimension
  CFuint _dim;
  
  /// stored the information about the setup that has been used (_set = 1, if setupOne has been used; _set = 2, if setupTwo has been used)
  CFuint _set;

  // ID of the starting cell
  CFuint _startCellID;
  
  CFuint _cellIdx;
  CFuint _FaceNb;
  
  /// ID of the exit face
  CFint _ExitFaceID;
  
  /// ID of the wall emitting face (used only if the emitting element is a wall and not a cell)
  CFuint _wallFaceID;
  
  /// Cell IDs map
  Common::CFMap<CFuint,CFuint> _CellIDmap;

}; // end of class ComputeNonLocalParameters

//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_ComputeNonLocalParameters_hh
