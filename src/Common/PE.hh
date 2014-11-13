// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Parallel_PE_hh
#define COOLFluiD_Parallel_PE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/PEInterface.hh"
#include "Common/WorkerStatus.hh"

#ifdef CF_HAVE_MPI
#  include "Common/MPI/PEInterfaceMPI.hh"
#else
#  include "Common/SERIAL/PEInterfaceSERIAL.hh"
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class controls the Parallel environment
class Common_API PE {
public:

    /// Initialise the PE
    static void InitPE (int * argc, char *** args);
    /// Checks if the PE is initialized
    static bool IsInitialised ();
    /// Return a reference to the current PE
    static PEInterface<>& GetPE ();
    /// Free the PE
    static void DonePE ();
  
#ifdef CF_HAVE_MPI
    /// Spawns processes on different host machines.
    /// @param count Number of processes to spawn.
    /// @param hosts Host list. Should be formatted as for -host option for mpirun
    /// (comma-separated items, without spaces - i.e. "host1,host2,host3"). May
    /// be null or empty (local host is used in these cases).
    /// @return Returns a MPI intercommunicator to allow the current process
    /// to communicate with created processes.
    /// @todo return a null comm or throw an exception on error (spawn failure,
    /// MPI env not initialized, count not valid...)
    static MPI_Comm spawn(unsigned int count, const char * hosts = CFNULL);
#endif
    
    /// Sets current process status.
    /// @param status New status
    static void setCurrentStatus(WorkerStatus::Type status);
    
    /// Gives the current process status.
    /// @return Returns the current process status
    static WorkerStatus::Type getCurrentStatus();
    
#ifdef CF_HAVE_MPI
     /// definition of class holding MPI group information
   class Group {
   public:
     
     /// default constructor
     Group() {}
     
     /// destructor
     ~Group() {MPI_Group_free(&group);}
     
     std::vector<int> globalRanks;
     std::vector<int> groupRanks; 
     MPI_Group group;
     MPI_Comm comm;
   };
   
   /// @return group corresponding to given global rank
   static CFint getGroupID(int rank) {return m_rank2Group.find(rank)->second;}
   
   /// @return the group data corresponding to given group ID
   static Group& getGroup(int groupID) 
   {
     cf_assert(static_cast<CFuint>(groupID) < m_groups.size());
     return *m_groups[groupID];
   }
   
   /// create MPI group 
   /// @param ranks          list of the ranks belonging to this group
   /// @param mapRank2Group  flag telling whether to build a reverse 
   ///                       rank-group mapping (each rank MUST be associated 
   ///                       to a unique group)
   static void createGroup(const std::vector<int>& ranks, 
			   const bool mapRank2Group); 
   
   /// clear the groups 
   /// @pre cannot be called from ~PE() because static PE object is destroyed 
   ///      after MPI_finalize(): it would try to double delete MPI_group's otherwise 
   static void clearGroups() {
     for (CFuint i = 0; i < m_groups.size(); ++i) { 
       delete m_groups[i]; 
     } 
   }
#endif
       
private:

    /// the current PE
    static PEInterface<> * m_curr_PE;

    /// Flag to keep track of Parallel Enviroment Initialized
    /// Note: cannot rely on m_curr_PE pointer because objects held by the enviroment
    ///       could not then check for initialization. An egg and chicken problem,
    ///       that appeared when using CFLog on the destructors of PEInterface related objects.
    static bool m_is_init;
    
    /// Path to the executable used to run the workers
    static char * m_command_workers;
    
    /// Current status. 
    /// Default value is @c #NOT_RUNNING.
    static WorkerStatus::Type m_current_status;
    
#ifdef CF_HAVE_MPI
    /// name of the object
    static std::map<CFint,CFint> m_rank2Group;
    
    /// list of groups
    static std::vector<Group*> m_groups; 
#endif
    
}; // end class PE

//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

#endif // COOLFluiD_Parallel_PE_hh
