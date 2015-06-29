// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PEInterface_hh
#define COOLFluiD_Common_PEInterface_hh

#include "Common/Group.hh"
#include "Common/PM.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Interface for the Parallel environment
/// @see PE
/// @author Dries Kimpe
template <typename PARMODEL = Common::PM_CUR>  class PEInterface;

/// Base class for the PEInterface
class Common_API PEInterfaceBase {
public:

  /// destructor
  ~PEInterfaceBase ();

  /// Returns the total number of execution contexts in
  /// the universum (not the big one of course, but the
  /// cluster-one ;-) )
  unsigned int GetProcessorCount () const;

  /// Return the ID of this processor (between 0 and GetProcessorCount)
  unsigned int GetRank () const;
  
  /// Set the barrier
  void setBarrier();

  /// Return true if this is a parallel simulation in some way
  /// (IMPORTANT: When COOLFluiD was compiled for multiple CPU's
  ///  - either MPI or SHM - but is started/called using only ONE
  ///  cpu, the run is NOT considered parallel!)
  bool IsParallel () const;

  /// returns the name of this model
  std::string GetName () const;

  /// return true if the model is multi-cpu capable
  bool IsParallelCapable () const;

  /// This function should be called periodically to help advance
  /// pending communication requests.
  void AdvanceCommunication ();
  
  /// @return group corresponding to given global rank
  std::string getGroupName(int rank) {return "Default";}
  
  /// @return the group data corresponding to given name
  Group& getGroup(const std::string name) {Group g; return g;}
  
  /// create MPI group 
  /// @param ranks          list of the ranks belonging to this group
  /// @param mapRank2Group  flag telling whether to build a reverse 
  ///                       rank-group mapping (each rank MUST be 
  ///                       associated to a unique group)
  void createGroup(const std::string nsp,
		   const std::string name,
		   const std::vector<int>& ranks, 
		   const bool mapRank2Group) {}
  
}; // end class PEInterface
      
//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_PEInterface_hh
