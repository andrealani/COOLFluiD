// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_MPIINITOBJECT_HH
#define PARALLEL_MPI_MPIINITOBJECT_HH

#include <mpi.h>

namespace COOLFluiD {
  namespace Common  {

//////////////////////////////////////////////////////////////////////////////

/// Interface to help other MPI dependent objects to initialize and cleanup
/// when MPI becomes available or is about to be unavailable.
class MPIInitObject
{
public:
    /// This function is called the moment MPI is fully initialized
    virtual void MPI_Init (MPI_Comm Communicator) = 0;

    /// This function is called *before* MPI deinitialization takes place.
    virtual void MPI_Done () = 0;

    virtual ~MPIInitObject ()
    {
    }
};

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif
