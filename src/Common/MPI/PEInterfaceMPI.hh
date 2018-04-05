// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PEInterfaceMPI_hh
#define COOLFluiD_Common_PEInterfaceMPI_hh

#include <list>
#include <map>

#include "Common/Group.hh"
#include "Common/PEInterface.hh"
#include "Common/NonCopyable.hh"
#include "Common/CFLog.hh"
#include "Common/MPI/MPIHelper.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD  {

  namespace Common {

    class MPIInitObject;
   
//////////////////////////////////////////////////////////////////////////////

/// Parallel Environment interface for MPI
/// Cannot be in the MPI namespace because it is a specialization of
/// a class in the Common namespace
template <>
class Common_API PEInterface<PM_MPI> : public PEInterfaceBase,
				       public Common::NonCopyable < PEInterface<PM_MPI> >
{
 public:
  
  /// Constructor (needs argc, args because of MPI)
  PEInterface (int * argc, char *** args);
  
  /// Destructor
  ~PEInterface ();
  
  /// Returns the total number of execution contexts in
  /// the universum (not the big one of course, but the
  /// cluster-one ;-) )
  unsigned int GetProcessorCount (const std::string nspaceName) const
  {
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::GetProcessorCount() => start\n");
    int nproc = 0;
    Common::CheckMPIStatus(MPI_Comm_size (GetCommunicator(nspaceName), &nproc));
    cf_assert(nproc > 0); 
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::GetProcessorCount() => end\n");
    return static_cast<unsigned int>(nproc);
  }
  
  /// Return the ID of this processor (between 0 and GetProcessorCount)
  unsigned int GetRank (const std::string nspaceName) const
  {
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::GetRank() => start\n");
    int rank = -1; 
    Common::CheckMPIStatus(MPI_Comm_rank (GetCommunicator(nspaceName), &rank));
    cf_assert(rank >= 0);
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::GetRank() => end\n");
    return static_cast<unsigned int>(rank);
  }
  
  /// advance communication
  /// This function should be called from time to time when there are
  /// overlapping
  /// communication requests (beginSync() is called but not endSync() )
  /// in orde to give the hardware/implementation a change to advance the
  /// communication.
  void AdvanceCommunication (const std::string nspaceName)
  {
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::AdvanceCommunication() => start\n");
    int Flag = 0;
    MPI_Comm comm = GetCommunicator(nspaceName);
    Common::CheckMPIStatus(MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, 
				       comm, &Flag, MPI_STATUS_IGNORE));
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::AdvanceCommunication() => end\n");
  }
  
  /// Return true if this is a parallel simulation in some way
  bool IsParallel() const 
  { 
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::IsParallel() => start\n"); 
    return (GetProcessorCount(std::string("Default")) > 1);
  }
  
  /// return true if the model is capable of multiple cpu's
  bool IsParallelCapable () const {return true;}
  
  /// return the name of the model
  std::string GetName () const {return "MPI";}
  
  /// init- en finalize functions
  /// The object will be freed by the PE (should
  /// become a ModulePointer)
  void RegisterInitObject (MPIInitObject * Obj);
  
  /// Set the barrier
  void setBarrier(const std::string nspaceName)  
  {
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::setBarrier() => start\n");
    MPI_Comm comm = GetCommunicator(nspaceName);  
    Common::CheckMPIStatus(MPI_Barrier(comm));
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::setBarrier() => end\n");
  }
  
  /// @return the communicator corresponding to the giving namespace name
  MPI_Comm GetCommunicator (const std::string nspaceName) const
  {
    CFLog(DEBUG_MIN, "PEInterface<PM_MPI>::GetCommunicator()\n");
    // if there are no groups, only the default exists 
    if (m_groups.size() == 0 || nspaceName == "Default" || m_groups.count(nspaceName) == 0) {
      return MPI_COMM_WORLD;
    }
    return m_groups.find(nspaceName)->second->comm;
  }
  
  /// @return group corresponding to given global rank
  std::string getGroupName(int rank) 
  {
    if (m_rank2Group.count(rank) == 0) {
      CFLog(ERROR, "PEInterfaceMPI::getGroupName() for rank " << rank << "not found!\n");
      abort();
    }
    return m_rank2Group.find(rank)->second;
  }
  
  /// @return the group data corresponding to given name
  Group& getGroup(const std::string& name)
  {
    if (!checkGroup(name)) {
      CFLog(ERROR, "PEInterfaceMPI::getGroup() => [" << name << "] does not exist!\n");
      cf_assert(checkGroup(name));
    }
    return *m_groups.find(name)->second;
  }
  
  /// @return a flag telling if the given name corresponds to an existing group
  bool checkGroup(const std::string& name)
  {
    return (m_groups.count(name) > 0);
  }
  
  /// create MPI group 
  /// @param nsp            namespace name
  /// @param name           group name (may be different from namespace)
  /// @param ranks          list of the ranks belonging to this group
  /// @param mapRank2Group  flag telling whether to build a reverse 
  ///                       rank-group mapping (each rank MUST be 
  ///                       associated to a unique group)
  void createGroup(const std::string nsp,
		   const std::string name,
		   const std::vector<int>& ranks, 
		   const bool mapRank2Group); 
  
  /// check if a given rank belongs to the specified group
  /// @param rank       rank to be checked for
  /// @param groupName  name of the group  which may (or not) contain the rank
  bool isRankInGroup(const int rank, const std::string groupName) 
  {
    if (groupName == "Default" || m_groups.count(groupName) == 0) return true;
    const Group& g = getGroup(groupName);
    return std::binary_search(g.globalRanks.begin(), g.globalRanks.end(), rank);
  }
  
  /// create the node to core mapping
  /// @param nsp            namespace name
  /// @param nodeIDs        array storing the node ID corresponding to the given nsp rank
  /// @param coreIDs        array storing the core ID corresponding to the given nsp rank
  /// @param uniqueNodeIDs  array of unique node IDs
  void createNodeToCoreMapping(const std::string nsp, 
			       std::vector<int>& nodeIDs,
			       std::vector<int>& coreIDs,
			       std::vector<int>& uniqueNodeIDs);
  
 private: // functions
  
  void CallInitFunctions ();
  void CallDoneFunctions ();
  void CallInitFunction (MPIInitObject * M) const;
  void CallDoneFunction (MPIInitObject * M) const;
    
  /// get the nodeID and coreID corresponding to the given rank 
  /// @author Wei Zhang, 2017-3-25, NSC
  void getID(MPI_Comm comm, const int size, const int rank, int& nodeID, int& coreID);
  
  /// get the nodeID and coreID corresponding to the given rank (new version, BG/Q compatible) 
  /// @author Wei Zhang, 2017-3-25, NSC
  void getIDNew(MPI_Comm comm, const int size, const int rank, int& nodeID, int& coreID);
  
  /// clear all groups 
  void clearGroups();
  
 private: // member data
  
  typedef std::list<MPIInitObject*> InitContainerType;
  InitContainerType InitList_;
  
  bool InitOK_;
  bool StopCalled_;
  
  /// name of the object
  std::map<int,std::string> m_rank2Group;
  
  /// list of groups
  std::map<std::string, Group*> m_groups; 
};
    
//////////////////////////////////////////////////////////////////////////////
    
  } // Common
  
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_PEInterfaceMPI_hh
