#ifndef COOLFluiD_Numerics_FluctSplit_PeriodicBCMPIeqX_hh
#define COOLFluiD_Numerics_FluctSplit_PeriodicBCMPIeqX_hh

//////////////////////////////////////////////////////////////////////////////
#include "mpi.h"
#include <map>
#include <iostream>
#include "Framework/DataSocketSink.hh"
#include "Framework/VectorialFunction.hh"

#include "Common/CFMap.hh"
#include "Common/PE.hh"

#include "FluctSplit/FluctuationSplitData.hh"



//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a periodic boundary condition for parralel computation
/// Connnect the nodes having the same X
/// Both boundaries of the periodic BC should be in the same TRS
/// @author Nadege Villedieu
class FluctSplit_API PeriodicBCMPIeqX : public FluctuationSplitCom {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  PeriodicBCMPIeqX(const std::string& name);

  /// Default destructor
  virtual ~PeriodicBCMPIeqX();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * synchronization of the processes
   */
  virtual void Barrier();

protected: // functions

  /// Execute on a set of dofs
  virtual void executeOnTrs();

  /**
   * Set the preProcesses connectivity between faces belonging to different process
   */
  virtual void preProcess();

protected: // data
 ///my process
  CFuint _my_nP;

  ///number processes
  CFuint _n_P;

  ///MPI_COMM_WORLD
  MPI_Comm _comm;

  /// map that maps global topological region face ID to a local one in the topological region set
  Common::CFMap<CFuint,CFuint> _globalToLocalTRSstateID;
  

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
   Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  ///number of faces per process
  std::vector<CFint> _countFpP;

/// map that maps a face with its corresponding periodic face using the index of the local process's trs
  Common::CFMap<CFuint,CFuint > _ConnectionStatePeriodic;

  /// map that maps a face with the number of the processor to which its corresponding periodic face belong
  Common::CFMap<CFuint,CFuint > _ConnectionProcessPeriodic;

  /// vector of send count for preProcess
  std::vector<int> _sendcounts2;

  /// vector of received count for preProcess
  std::vector<int> _recvcounts2;

  /// vector of received displacements for preProcess
  std::vector<int> _rdispls2;

  /// vector of sed displacements for preProcess
  std::vector<int> _sdispls2;

  /// vector of send count for preProcess
  std::vector<int> _sendcounts3;

  /// vector of received count for preProcess
  std::vector<int> _recvcounts3;

  /// vector of send count for preProcess
  std::vector<int> _sendcounts4;

  /// vector of received count for preProcess
  std::vector<int> _recvcounts4;

  /// vector of received displacements for preProcess
  std::vector<int> _rdispls3;

  /// vector of sed displacements for preProcess
  std::vector<int> _sdispls3;

  ///number of equation
  CFuint _nE;

  /// vector of received buffer for pre Process
  std::vector<double> _rbuf2;

  /// vector of displacements between the arrays of boundary state of the faces belonging to a single process
  std::vector<CFuint> _RecvDis2;

  /// vector of received buffer for pre Process
  std::vector<double> _rbuf3;

  /// vector of displacements between the arrays of boundary state of the faces belonging to a single process
  std::vector<CFuint> _RecvDis3;

  /// vector of received buffer for pre Process
  std::vector<double> _rbuf4;

  /// vector of displacements between the arrays of boundary state of the faces belonging to a single process
  std::vector<CFuint> _RecvDis4;

  /// vector of boundary state of all faces belonging to a single process
  std::vector<CFreal> _Boundaryrhs;

  /// vector of boundary state of all faces belonging to a single process
  std::vector<CFreal> _Boundaryupcoef;

  /// vector of boundary state of all faces belonging to a single process
  std::vector<CFreal> _Boundaryupcoef2;

  /// the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// number of nodes
  CFuint m_nb_st_trs;

 /// size of received buffer
  CFuint _LastDisplacement;

 /// size of received buffer
  CFuint _LastDisplacement3;

  // Framework::DataHandle< Framework::State*, Framework::GLOBAL > m_states;
  /// the socket to the data handle of the isUpdated flags
  // Framework::DataSocketSink<bool> socket_isUpdated;

  /// name of the trs to be coupled
  std::string m_coupled_trs;

  /// vector holding the indexes of the states on the coupled TRS that match the
  /// applied TRS
  // std::vector < CFuint > m_match_states_idx;

  /// transformed coordinate
  // RealVector m_tcoord;



  /// difference between transformed and coupled state coordinates
  // RealVector m_delta;

  /// temporary residual
  // RealVector m_tmp_rhs;

  /// threshold of distance that will consider two states matching after the transformation of coordinates
  CFreal m_threshold;

  /// socket for isUpdated
  Framework::DataSocketSink< bool> socket_isUpdated;

  /// a vector of string to hold the functions for transformation of coordinates
  //std::vector<std::string> m_transform_funcs;

  /// a vector of string to hold the variables
  //std::vector<std::string> m_vars;

  /// parser for the functions for transformation of coordinates
  //Framework::VectorialFunction m_vFunction;

}; // end of class PeriodicBCMPIeqX

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PeriodicBCMPIeqX_hh
