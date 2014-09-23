// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/MPI/MPIHelper.hh"
#include "Common/MPI/MPIException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

void DoMPIError (int status)
{
    char Buf[MPI_MAX_ERROR_STRING+1];
    int BufSize = sizeof (Buf)-1;
    MPI_Error_string (status, Buf, &BufSize);
    std::string S ("MPI Error: ");
    S+=Buf;
    S+='\n';
    throw MPIException (FromHere(),S);
}

//////////////////////////////////////////////////////////////////////////////

    }
}
