#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCNeumannFromFile_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCNeumannFromFile_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Common/CFMap.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/StateInterpolator.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"

#include "Common/ConnectivityTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace Common {
    template <class KEY, class VALUE> class LookUpTable;
  }
    
    namespace Framework {
    class Node;
    class State;
    class TopologicalRegionSet;
  }
    
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Neumann boundary condition with data from an input file
 *
 * @author Ray Vandenhoeck
 */
class BCNeumannFromFile : public BCStateComputer {

public:  // methods
    
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCNeumannFromFile(const std::string& name);

  /// Destructor
  ~BCNeumannFromFile();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNeumannFromFile";
  }

  /// Setup private data
  void setup();

  /// Unsetup private data
  void unsetup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
   void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                              std::vector< std::vector< RealVector* > >& ghostGrads,
                              const std::vector< RealVector >& normals,
                              const std::vector< RealVector >& coords);
   
   
   protected: // helper function
  
  /**
   * Local nested class for grouping surface data
   */
  class SurfaceData {
  public:
    RealMatrix xyz;
    RealVector Tw;
  };
  
  /**
   * Local nested class for grouping data related to the closet point
   */
  class ClosestPointData {
  public:
    std::vector<CFint> surfaceIDs;
    std::vector<CFint> pointsIDs;
    RealVector r;
    
    void reset() 
    {
      surfaceIDs.assign(surfaceIDs.size(), -1); 
      pointsIDs.assign(pointsIDs.size(), -1);
      r = MathTools::MathConsts::CFrealMax();
    }
    
    void regressionFromTo(CFuint start, CFuint end)
    {
      cf_assert(end < surfaceIDs.size());
      cf_assert(start < surfaceIDs.size());
      surfaceIDs[end] = surfaceIDs[start]; 

      cf_assert(end < pointsIDs.size());
      cf_assert(start < pointsIDs.size());
      pointsIDs[end] = pointsIDs[start];

      cf_assert(end < r.size());
      cf_assert(start < r.size());
      r[end] = r[start];
    }
  };
  
  /**
   * Read the surface data
   */
  void readSurfaceData(std::vector<SurfaceData*>& surfaces);
  
  /**
   * Read line data at z=0 from given surface
   */
  SurfaceData* extractLineData(SurfaceData* surface);
  
  protected: // data
  
  /// ID corresponding to the z coordinate (z=0) for which plane is extracted
  CFint m_extractCoordZID;
  
  // name of the TRSs on which values must be prescribed
  std::vector<std::string> m_trsName;
  
  /// name of the file where the temperature distribution is provided
  std::string m_fileNameTw;
  
  /// rotation angle
  CFreal m_angle;
  
  /// ID  of the spatial coordinates lying in the rotation plane
  std::vector<CFuint> m_xvec;
  
  /// add minus sign to read data
  bool m_addMinus;
  
  /// number of closest points for interpolation stencil
  CFuint m_nbClosestPoints;
  
  /// IDs corresponding to the x,y coordinate (z=0) for which plane is extracted
  std::vector<CFint> m_extractCoordXYID;
  
    // flx pnt mapped coordinates
  std::vector< RealVector > m_flxPntsLocalCoords;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  
  /**
   * Compute the periodic state given the face. 
   * It is assumed that MPI-communication has occured in preProcess() first
   */
  //Framework::State* computePeriodicState(Framework::GeometricEntity *const face);
  //till here
  /// TRS in which this BC is defined
  Common::SafePtr< Framework::TopologicalRegionSet > m_thisTRS;

  /// local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;

  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;  
  
  /// number of dimensions in the physical model
  CFuint m_dim;

  /// variable for current face orientation
  CFuint m_orient;

  /// face builder
  //Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
  //Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> _faceBuilder;
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > _faceBuilder;

    
  /// map that maps global topological region face ID to a local one in the topological region set
  Common::CFMap<CFuint,CFuint> m_globalToLocalTRSFaceID;

  /**
   * Definition of a structure containing face-centred coordinates and faceID's of a face
   */
  class FlxPntStruct
  {
    public:
      RealVector _centreCoordinates;
      int     _globalFaceID;
      int     _localFaceID;
      int     _localFluxID;
      int     _orient;
      CFreal     _Tw;

    public:
      void setOrient(CFuint orient){ _orient=orient;}
      CFuint getOrient() const {return _orient;}
      void setGlobalFaceID(CFuint globalFaceID) { _globalFaceID = globalFaceID;}
      CFuint getGlobalFaceID() const {return _globalFaceID;}
      void setLocalFaceID(CFuint localFaceID) { _localFaceID = localFaceID;}
      CFuint getLocalFaceID() const {return _localFaceID;}
      void setLocalFluxID(CFuint localFluxID) { _localFluxID = localFluxID;}
      CFuint getLocalFluxID() const {return _localFluxID;}
      void setCentreCoordinates(const RealVector& centreCoordinates) { _centreCoordinates.resize(centreCoordinates.size()); _centreCoordinates = centreCoordinates;}
      const RealVector& getCentreCoordinates() const {return _centreCoordinates;}
      void setTw(CFreal Tw){ _Tw=Tw;}
      CFreal getTw() const {return _Tw;}
  };
  
  // boundary flux points
  RealVector m_flxPntTws;
  
  //No. of faces in the TRS
  CFuint m_nbGeoEnts;
  
  // number of equations
  CFuint m_nbEqs;

  protected: // method
  
  void createFaceOrientationStartIndexes();

}; // class BCNeumann

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNeumannFromFile_hh

