// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_MPI_MPIError_hh
#define COOLFluiD_Common_MPI_MPIError_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>
#include <string>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// The following class reports errors from MPI calls
/// @author Andrea Lani
class MPIError {
public:
  
  /// get instance of his singleton
  static MPIError& getInstance() 
  {
    static MPIError mErr;
    return mErr;
  }
  
  /// initialize the error handler
  void init(MPI_Comm comm, int rank) 
  {
    m_comm = comm;
    m_rank = rank;
    m_errorCode = MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
  }
  
  /// check the sanity of the MPI call
  void check(const std::string& mpiFun, const std::string& fromHere, int errorCode)
  {
    using namespace std;
    
    m_errorCode = errorCode;
    int msgLength = 0;
    if (m_errorCode != MPI_SUCCESS) {
      MPI_Error_string(m_errorCode, m_errorMsg, &msgLength);
      cout << "P"<< m_rank << "=> Error in call to " << mpiFun << ": " << string(m_errorMsg) << "\n";
      cout << "P"<< m_rank << "=> Exiting from " << fromHere << "\n";
      MPI_Abort(m_comm, -1);
    }
  }
  
private:
  
  /// communicator
  MPI_Comm m_comm;
   
  /// rank of the current processor
  int m_rank;
  
  /// error code
  int m_errorCode;
   
  /// error message
  char m_errorMsg[MPI_MAX_ERROR_STRING];
  
};
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_MPI_MPIStructDef_hh
