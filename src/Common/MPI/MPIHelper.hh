// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_MPIHELPER_HH
#define PARALLEL_MPI_MPIHELPER_HH

#include <mpi.h>

#include "Common/CommonAPI.hh"

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Act on an error
Common_API
void DoMPIError (int status);

/// Lightweight & fast check to see if there is an error
/// Calls DoMPIError if there is one
/// Once proper support for throwing errors is confirmed, this check go.
Common_API
inline void CheckMPIStatus (int status)
{
  if ( status != MPI_SUCCESS ) DoMPIError (status);
}

//////////////////////////////////////////////////////////////////////////////

   }
}

#endif
