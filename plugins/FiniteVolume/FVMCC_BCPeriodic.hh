#ifndef COOLFluiD_Numerics_FiniteVolume_BCPeriodic_hh
#define COOLFluiD_Numerics_FiniteVolume_BCPeriodic_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/CFMap.hh"
#include "Common/PE.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a parallel periodic boundary condition command for
   * two topological regions in 1D, 2D and 3D.
   * 
   * Topological regions have to be combined in 1 TRS
   * A Translation Vector must be given between the 2 regions
   *
   * @author Willem Deconinck
   */
   
class BCPeriodic : public FVMCC_BC {
public:

   /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  BCPeriodic(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BCPeriodic();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Transfer the gradients data
   */
  virtual void transferGradientsData();
  
  /**
   * Set the preProcesses connectivity between faces belonging to different processes
   */
  virtual void preProcess();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
   
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /// Compute the periodic cell gradient and corresponding limiter
  virtual void computePeriodicGradient(Framework::GeometricEntity *const face, 
				       std::vector<CFreal*>& gradientLimiter,
				       CFreal*& cellLimiter);
  
  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "BCPeriodic";
  }
  
protected:
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder; 
  
  /// setup MPI parameters
  virtual void setupMPI();
  
  /**
   * Compute the periodic state given the face. 
   * It is assumed that MPI-communication has occured in preProcess() first
   */
  Framework::State* computePeriodicState(Framework::GeometricEntity *const face);
  
  /**
   * Return if the face is an outlet face. (outlet = true, inlet = false)
   */
  bool isOutletFace(Framework::GeometricEntity *const face)
  {
    const CFuint faceGlobalID = face->getID();
    const CFuint faceLocalID = _globalToLocalTRSFaceID.find(faceGlobalID);
    return _localWestEastMap[faceLocalID];
  }
  
  bool isOutletFace(const CFuint faceLocalID)
  {
    return _localWestEastMap[faceLocalID];
  }
  
  /**
   * Return if the face is an inlet face. (inlet = true, outlet = false)
   */
  bool isInletFace(Framework::GeometricEntity *const face)
  {
    const CFuint faceGlobalID = face->getID();
    const CFuint faceLocalID = _globalToLocalTRSFaceID.find(faceGlobalID);
    return !_localWestEastMap[faceLocalID];
  }
  
  bool isInletFace(const CFuint faceLocalID)
  {
    return !_localWestEastMap[faceLocalID];
  }
 
  /**
   * Number of periodic faces
   */
   CFuint getNbPeriodicFaces()
   {
     return _nbTrsFaces;
   }
  
  /// transfer the gradient data
  void transferArray(const Framework::DataHandle<CFreal>& gradient,
		     std::vector<CFreal>& sendbuf, 
		     std::vector<CFreal>& recvbuf);
  
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
  class FaceStruct
  {
    public:
      RealVector _centreCoordinates;
      int     _globalFaceID;
      int     _localFaceID;
      
    public:
      void setGlobalFaceID(CFuint globalFaceID) { _globalFaceID = globalFaceID;}
      CFuint getGlobalFaceID() const {return _globalFaceID;}
      void setLocalFaceID(CFuint localFaceID) { _localFaceID = localFaceID;}
      CFuint getLocalFaceID() const {return _localFaceID;}
      void setCentreCoordinates(const RealVector& centreCoordinates) { _centreCoordinates.resize(centreCoordinates.size()); _centreCoordinates = centreCoordinates;}
      const RealVector& getCentreCoordinates() const {return _centreCoordinates;}
  };
  
  /**
   * FaceStruct with MPI functionality
   */ 
  class FaceMPIStruct : public FaceStruct {
  private:
    Common::MPIStruct ms;
    int ln[3];
    MPI_Comm _comm;
    
  public:
    FaceMPIStruct() : FaceStruct() { 
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
    
    void copy(const FaceStruct& face) {
      setLocalFaceID(face.getLocalFaceID());
      setGlobalFaceID(face.getGlobalFaceID());
      setCentreCoordinates(face.getCentreCoordinates());
    }
  };
  
  /**
   * Functor used in the std::find_if algorithm to match 2 periodic vectors
   */
  struct findTranslated 
  {
    RealVector _Vt;
    const CFreal _threshold;
    
    findTranslated(const FaceStruct& firstVector, const RealVector& translationVector, const CFreal& threshold) : _threshold(threshold)
    {
      _Vt.resize(translationVector.size());
      _Vt = firstVector.getCentreCoordinates() + translationVector; 
    }
    
    bool operator()(const FaceStruct& iter)
    {
      bool isMatch = true;
      const RealVector& iterVec = iter.getCentreCoordinates();
      const CFuint dim = iterVec.size();
      //CFLog(INFO,"Translated Coords"<<_Vt<<"\tCompared coords = "<<iterVec<<"\n");
      for(CFuint i=0; i<dim; ++i) {
        isMatch = isMatch && (std::abs(iterVec[i]-_Vt[i]) < _threshold);
      }
      return isMatch;
    }
  };

protected: //data
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// socket for uZ values
  Framework::DataSocketSink<CFreal> socket_uZ;
  
  /// socket for limiter
  Framework::DataSocketSink<CFreal> socket_limiter;
  
  /// flag telling if socket uX exists
  bool _uxExists;
  
  /// flag telling if socket uY exists
  bool _uyExists;
  
  /// flag telling if socket uZ exists
  bool _uzExists;
  
  /// Translation Vector between planes
  std::vector<CFreal> _translationVector;
  
  /// map that maps global topological region face ID to a local one in the topological region set
  Common::CFMap<CFuint,CFuint> _globalToLocalTRSFaceID;
  
  /// global west or east flag  (west = 0, east = 1)
  Common::CFMap<CFuint,bool> _localWestEastMap;
  
  /// threshold
  CFreal _threshold;
 
  /// number of equation
  CFuint _nE;

  /// number of faces
  CFuint _nbTrsFaces;
  
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
  
}; // end of class BCPeriodic

//////////////////////////////////////////////////////////////////////////////


 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BCPeriodic_hh
