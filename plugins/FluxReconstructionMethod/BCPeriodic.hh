#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCPeriodic_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCPeriodic_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "Common/MPI/MPIStructDef.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents periodic boundary conditions
 *
 * @author Ray Vandenhoeck
 * @author Isaac Alonso Asensio
 */
class BCPeriodic : public BCStateComputer {

public:  // methods

  /// Constructor
  BCPeriodic(const std::string& name);

  /// Destructor
  ~BCPeriodic();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCPeriodic";
  }

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);


  /// Set up private data
  void setup();
  
  /// Unset up private data
  void unsetup();
  
  //Added from FV
   /**
   * Set the preProcesses connectivity between faces belonging to different processes
   */
  virtual void preProcess();
  //till here
  

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

  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  void createFaceOrientationStartIndexes();

protected:
    //Added from FV
    
    // flx pnt mapped coordinates
  std::vector< RealVector > m_flxPntsLocalCoords;
  
  // All flux points of a cell
  Common::SafePtr<std::vector< RealVector > > m_allCellFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;
  
  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< CFint > > m_faceMappedCoordDir;
  
  /// vector for the face jacobian vectors
  std::vector< RealVector > m_faceJacobVecs;
  
 
  
   /// face Jacobian vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;
  
  /// face Jacobian vector sizes
  std::vector< CFreal > m_faceJacobVecSizeFlxPnts;
  
   /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;

    /// setup MPI parameters
  virtual void setupMPI();
  
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

  /// nbr of sol pnts on which a flx pnt is dependent
  CFuint m_nbrSolDep;

  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;

  /// Boundary face at the other side of the domain
  Framework::GeometricEntity* m_otherFace;

  /// variable for current cell
  Framework::GeometricEntity* m_intCell;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// the gradients in the neighbouring cell
  std::vector< std::vector< RealVector >* > m_cellGrads;

  /// Maps each flux point to their correspoding periodic face
  std::vector< CFuint > _faceConnectivityMap;

  /// Maps each flux point to their corresponding periodic local flux point
  std::vector< CFuint > _fluxPointConnectivityMap;

  /// Collect the orientation of the face at the other side of the boundary
  std::vector< CFuint > _orientMap;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDep;

  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;

  /// face builder
  //Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
  //Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> _faceBuilder;
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > _faceBuilder;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
   /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;

  //Framework::DataSocketSink<CFreal> socket_normals;
  //Framework::DataSocketSink<CFreal> socket_states;
  //Framework::DataSocketSink<CFreal> socket_gstates;
  //Framework::DataSocketSink<CFreal> socket_nodes;
    
  /// map that maps global topological region face ID to a local one in the topological region set
  Common::CFMap<CFuint,CFuint> _globalToLocalTRSFaceID;


  /// global west or east flag  (west = 0, east = 1)
  Common::CFMap<CFuint,bool> _localWestEastMap;
  
  /// threshold
  CFreal _threshold;
        
  /**
   * Definition of a pair containing the processor number and the faceID
   */
  class PairStruct
  {
    private:
      CFuint _process;
      CFuint _faceID;
    public:
      PairStruct(const CFuint& process, const CFuint& faceID) : _process(process), _faceID(faceID) {}
      CFuint getFaceID() {return _faceID;}
      CFuint getProcess() {return _process;}
  };

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
  };
  
  /**
   * FaceStruct with MPI functionality
   */ 
  class FaceMPIStruct : public FlxPntStruct {
  private:
    Common::MPIStruct ms;
    int ln[3];
    MPI_Comm _comm;
    
  public:
    FaceMPIStruct() : FlxPntStruct() { 
      CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
      ln[0] = 1;     // first block is localFaceID --> 1 element
      ln[1] = 1;     // second block is globalFaceID --> 1 element
      ln[2] = dim;   // third block is Face coordinates --> "dim" elementsbuild(); 
      _centreCoordinates.resize(dim);
      
      const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
      _comm = Common::PE::GetPE().GetCommunicator(nsp);
      Common::MPIStructDef::buildMPIStruct(&_localFaceID, &_globalFaceID, &_centreCoordinates[0], ln, ms);
    }
    
    void broadcast(const CFuint& iP) {
      MPI_Bcast(ms.start, 1, ms.type, iP, _comm);
    }
    
    void copy(const FlxPntStruct& face) {
      setLocalFaceID(face.getLocalFaceID());
      setGlobalFaceID(face.getGlobalFaceID());
      setCentreCoordinates(face.getCentreCoordinates());
      setLocalFluxID(face.getLocalFaceID());
    }
  };


  /**
 *    * Functor used in the std::find_if algorithm to match 2 periodic vectors
 *       */
  struct findTranslated
  {
    RealVector _Vt;
    const CFreal _threshold;

    findTranslated(const FlxPntStruct& firstVector, const RealVector& translationVector, const CFreal& threshold) : _threshold(threshold)
    {
      _Vt.resize(translationVector.size());
      _Vt = firstVector.getCentreCoordinates() + translationVector;
    }

    bool operator()(const FlxPntStruct& iter)
    {
      bool isMatch = true;
      const RealVector& iterVec = iter.getCentreCoordinates();
      const CFuint dim = iterVec.size();
      //CFLog(INFO,"Translated Coords "<<_Vt<<"\tCompared coords = "<<iterVec<<"\n");
      for(CFuint i=0; i<dim; ++i) {
        isMatch = isMatch && (std::abs(iterVec[i]-_Vt[i]) < _threshold);
      }
      return isMatch;
    }
  };

  /// Translation Vector between planes
  std::vector<CFreal> _translationVector;
  //Added from FV
protected: //data
    
    //No. of equations
    CFuint NbEqs;
    
    
    //No. of faces in the TRS
    CFuint nbGeoEnts;
    
    
    /// maps local face to periodic process and periodic face (can be more than 1)
  std::vector< std::vector< PairStruct > > _localConnectivityMap;
    
    /**
   * MPI Related parameters
   */
   
  /// MPI_COMM_WORLD
  MPI_Comm _comm;

  /// number of processes
  CFuint _nbProcesses;
  
  /// rank of this process
  CFuint _rank;
  
  /// number of faces per process
  std::vector<CFuint> _nbFacesPerProcess;

  /// A temporary container for a state
  Framework::State* _periodicState;

  /// SendReceive parameters
  std::vector<int> _sendcounts;
  std::vector<int> _recvcounts;
  std::vector<int> _recvdispls;
  std::vector<int> _senddispls;
  std::vector<double> _recvbuf;
  std::vector<double> _sendbuf;

  std::vector<double> _recvbufLimiter;
  std::vector<double> _sendbufLimiter;
  std::vector<std::vector<double> > _recvbufGrad;
  std::vector<std::vector<double> > _sendbufGrad;
    
  CFuint _LastDisplacement;
  //till here
  
}; // class BCPeriodic
    
    //////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCPeriodic_hh

