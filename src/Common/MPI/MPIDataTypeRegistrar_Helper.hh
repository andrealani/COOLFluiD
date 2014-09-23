// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MPIDATATYPEREGISTRAR_HELPER_HH
#define MPIDATATYPEREGISTRAR_HELPER_HH

#include "Common/MPI/MPIDataType.hh"
#include <mpi.h>

namespace COOLFluiD {
    namespace Common  {
//////////////////////////////////////////////////////////////////////////////

class MPIDataTypeRegistrar_Helper : public MPIDataType
{
protected:
    MPI_Datatype TheType;

public:
    MPIDataTypeRegistrar_Helper ();
    virtual ~MPIDataTypeRegistrar_Helper ();


    void DoRegister ();

    MPI_Datatype GetType () const
    { return TheType; }

};

//////////////////////////////////////////////////////////////////////////////

   }
}

#endif
