// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CommPattern.hh"

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

CommPattern::CommPattern (MPI_Comm C)
    : Communicator_ (C)
{
    MPI_Comm_size (Communicator_, &Size_);
    MPI_Comm_rank (Communicator_, &Rank_);
}

MPI_Comm CommPattern::GetCommunicator () const
{
    return Communicator_;
}


//////////////////////////////////////////////////////////////////////////////

    }
}
