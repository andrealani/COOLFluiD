#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCInletHelioUnsteadyMHD_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCInletHelioUnsteadyMHD_hh

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

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an Inlet boundary condition
 * for the 3D MHD equations.
 *
 * @author Rayan Dhib
 */
class BCInletHelioUnsteadyMHD : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCInletHelioUnsteadyMHD(const std::string& name);

  /// Destructor
  ~BCInletHelioUnsteadyMHD();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCInletHelioUnsteadyMHD";
  }

  /// Set up private data and data
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
             
   /**
   * Set the B boundary values
   */
   virtual void preProcess();


protected: // helper function
  
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
      CFreal     _rho;
      CFreal     _u;
      CFreal     _v;
      CFreal     _w;
      CFreal     _Bx;
      CFreal     _By;
      CFreal     _Bz;
      CFreal     _p;

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
      void setRhow(CFreal rho){ _rho=rho;}
      void setUw(CFreal u){ _u=u;}
      void setVw(CFreal v){ _v=v;}
      void setWw(CFreal w){ _w=w;}
      void setBxw(CFreal Bx){ _Bx=Bx;}
      void setByw(CFreal By){ _By=By;}
      void setBzw(CFreal Bz){ _Bz=Bz;}
      void setPw(CFreal p){ _p=p;}
      CFreal getRhow() const {return _rho;}
      CFreal getUw() const {return _u;}
      CFreal getVw() const {return _v;}
      CFreal getWw() const {return _w;}
      CFreal getBxw() const {return _Bx;}
      CFreal getByw() const {return _By;}
      CFreal getBzw() const {return _Bz;}
      CFreal getPw() const {return _p;}
  };
  
  /**
   * Local nested class for grouping surface data
   */
  class SurfaceData {
  public:
    RealMatrix xyz;
    RealVector rho;
    RealVector u;
    RealVector v;
    RealVector w;
    RealVector Bx;
    RealVector By;
    RealVector Bz;
    RealVector p;
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
  void readSurfaceData(std::vector<SurfaceData*>& surfaces, const std::string& fileName);
  
  /**
   * Read line data at z=0 from given surface
   */
  SurfaceData* extractLineData(SurfaceData* surface);

  void createFaceOrientationStartIndexes();
  
  /**
   * Interpolate surface data between files at current time
   */
  void interpolateSurfaceDataInTime();
  
  /**
   * Dynamic time interpolation for specific time (called from computeGhostStates)
   */
  void interpolateSurfaceDataAtTime(CFreal targetTime);
  
  /**
   * Update boundary flux point values after time interpolation
   */
  void updateBoundaryFluxPointValues();
  
  /**
   * Extrapolate spatial data from current time-interpolated surface
   */
  void extrapolateSpatialData();
  
  /**
   * Perform initial spatial interpolation for single file case
   */
  void performInitialSpatialInterpolation(const std::vector<SurfaceData*>& surfaces, 
                                          const std::vector<FlxPntStruct>& bndFlxPnts);
  
  /**
   * Generate file names and times automatically based on pattern and time range
   */
  void generateFileListFromPattern();
  

protected: // data

  /// map TRS name -> initial solution array that will be used as BC value
  Common::CFMap<std::string, RealVector*> m_initialSolutionMap;
  
  /// array specifying IDs of initial solution components that will be used as BC value
  std::vector<CFuint> m_initialSolutionIDs;

  /// variable for current internal cell
  Framework::GeometricEntity* m_intCell;

  /// current face
  Framework::GeometricEntity* m_currFace;

  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;
    
  /// ID corresponding to the z coordinate (z=0) for which plane is extracted
  CFint m_extractCoordZID;
  
  // name of the TRSs on which values must be prescribed
  std::vector<std::string> m_trsName;
  
  /// name of the file where the wall values is provided
  std::string m_fileName;
  
  /// rotation angle
  CFreal m_angle;
  
  /// ID  of the spatial coordinates lying in the rotation plane
  std::vector<CFuint> m_xvec;
  
  /// number of closest points for interpolation stencil
  CFuint m_nbClosestPoints;
  
  /// IDs corresponding to the x,y coordinate (z=0) for which plane is extracted
  std::vector<CFint> m_extractCoordXYID;
  
  /// flx pnt mapped coordinates
  std::vector< RealVector > m_flxPntsLocalCoords;
  
  /// TRS in which this BC is defined
  Common::SafePtr< Framework::TopologicalRegionSet > m_thisTRS;

  /// local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;

  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;
  
  /// local coordinates of the flux points on one face per face type
  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerType;
  
  /// number of dimensions in the physical model
  CFuint m_dim;

  /// variable for current face orientation
  CFuint m_orient;

  /// face builder for boundary conditions
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > _faceBuilder;
    
  /// map that maps global topological region face ID to a local one in the topological region set
  Common::CFMap<CFuint,CFuint> m_globalToLocalTRSFaceID;

  // boundary flux points
  RealVector m_flxPntRhos;
  RealVector m_flxPntUs;
  RealVector m_flxPntVs;
  RealVector m_flxPntWs;
  RealVector m_flxPntBxs;
  RealVector m_flxPntBys;
  RealVector m_flxPntBzs;
  RealVector m_flxPntPs;
  
  //No. of faces in the TRS
  CFuint m_nbGeoEnts;
  
  // number of equations
  CFuint m_nbEqs;

  /// physical model
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> m_varSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;

  /// Max number of flux pnts on a face
  CFuint m_nbrFaceFlxPntsMax;
  
  /// inner states
  std::vector< RealVector > m_tempStates;

    /// file names for time interpolation
  std::vector<std::string> m_fileNameTw;
  
  /// times corresponding to each file
  std::vector<CFreal> m_fileNameTime;
  
  /// all surface data for time interpolation
  std::vector<std::vector<SurfaceData*> > m_allSurfaces;
  
  /// current interpolated surface data
  std::vector<SurfaceData*> m_surfaceAtTime;
  
  /// flag to check if time interpolation is needed
  bool m_useTimeInterpolation;
  
  // Automatic file generation parameters
  /// pattern for automatic file generation
  std::string m_filePattern;
  
  /// enable automatic file generation
  bool m_autoGenerate;
  
  /// starting date-time string
  std::string m_startDateTime;
  
  /// ending date-time string  
  std::string m_endDateTime;
  
  /// time step between files in non-dimensional units
  CFreal m_dataTimeStep;
  
  /// starting simulation time (non-dimensional)
  CFreal m_startTime;
  
  /// ending simulation time (non-dimensional)
  CFreal m_endTime;
  
  /// last time when interpolation was performed (for efficiency)
  CFreal m_lastInterpolationTime;

}; // class BCInletHelioUnsteadyMHD

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCInletHelioUnsteadyMHD_hh

