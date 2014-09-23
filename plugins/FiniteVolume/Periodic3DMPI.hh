#ifndef COOLFluiD_Numerics_FiniteVolume_PeriodicX3DMPI_hh
#define COOLFluiD_Numerics_FiniteVolume_PeriodicX3DMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/CFMap.hh"
#include "mpi.h"
#include "Common/PE.hh"
#include <map>


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a periodic boundary condition command for two
   * topological regions in 3D parallel to x-axis for CellCenterFVM schemes
   * for serial simulations
   *
   * @author Alessandro Sanna
   *
   */
class Periodic3DMPI : public FVMCC_BC {
public:

   /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  Periodic3DMPI(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~Periodic3DMPI();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Set the preProcesses connectivity between faces belonging to different process
   */
  virtual void preProcess();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

  /**
   * synchronization of the processes
   */
  virtual void Barrier();

private: //data

  /// map that maps global topological region face ID to a local one in the topological region set
  Common::CFMap<CFuint,CFuint> _globalToLocalTRSFaceID;
  
  /// map that maps a face with its corresponding periodic face using the index of the local process's trs
  Common::CFMap<CFuint,CFuint > _ConnectionFacePeriodic;

  /// map that maps a face with the number of the processor to which its corresponding periodic face belong
  Common::CFMap<CFuint,CFuint > _ConnectionProcessPeriodic;
  
  /// vector of displacements between the arrays of boundary state of the faces belonging to a single process
  std::vector<CFuint> _RecvDis2;

  /// vector of send count for preProcess
  std::vector<int> _sendcounts2;

  /// vector of received count for preProcess
  std::vector<int> _recvcounts2;

  /// vector of received displacements for preProcess
  std::vector<int> _rdispls2;

  /// vector of sed displacements for preProcess
  std::vector<int> _sdispls2;

  /// size of received buffer
  CFuint _LastDisplacement;

  /// vector of received buffer for pre Process
  std::vector<double> _rbuf2;

  /// vector of boundary state of all faces belonging to a single process
  std::vector<CFreal> _BoundaryState;

  ///number of faces per process
  std::vector<CFint> _countNWf;

  ///number of faces per process
  std::vector<CFint> _countNEf;

  ///number of faces per process
  std::vector<CFint> _countFpP;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder; 
 
  ///number of equation
  CFuint _nE;

  ///number of faces
  CFuint _nbTrsFaces;

  ///my process
  CFuint _my_nP;

  ///number processes
  CFuint _n_P;

  ///MPI_COMM_WORLD
  MPI_Comm _comm;
  
  /// threshold
  CFreal _threshold;
  
}; // end of class Periodic3DMPI

//////////////////////////////////////////////////////////////////////////////


 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Periodic3DMPI_hh
