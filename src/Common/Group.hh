// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Group_hh
#define COOLFluiD_Common_Group_hh

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
#include <mpi.h>
#include "MPI/MPIHelper.hh"
#endif

#include <vector>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Class holding MPI group information
/// @author Andrea Lani
class Group {
public:

  /// default constructor
  Group() {}
  
  /// destructor
  ~Group() 
  {
#ifdef CF_HAVE_MPI
    CheckMPIStatus(MPI_Group_free(&group));
#endif
  }
  
  // global MPI ranks associated to this group
  std::vector<int> globalRanks;
  
  // MPI ranks within the group
  std::vector<int> groupRanks; 
  
#ifdef CF_HAVE_MPI
  // group object 
  MPI_Group group;
  
  // group communicator 
  MPI_Comm comm;
#endif
  
}; // end class Group
      
//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_PE_hh
