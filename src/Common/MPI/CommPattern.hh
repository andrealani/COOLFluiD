// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_COMMPATTERN_HH
#define PARALLEL_MPI_COMMPATTERN_HH

#include <mpi.h>

#include "Common/COOLFluiD.hh"

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// base class for communication patterns
///   stores communicator and constants
class CommPattern
{
public:

    // value type: first=receive from, second=send to
    typedef std::pair<int,int> value_type;


    /// Constructor
    CommPattern (MPI_Comm C = MPI_COMM_SELF);

    /// Return the communicator to be used in communications
    /// (could be different from the one the pattern was initialised with
    MPI_Comm GetCommunicator () const;


    /// Return the communicator size
    int GetCommSize () const
    {
return Size_;
    }

    /// Return our rank
    int GetRank () const
    {
return Rank_;
    }

protected:
    MPI_Comm Communicator_;

    int Rank_;
    int Size_;
};

template <typename ITERATOR>
class CommPatternBase : public CommPattern
{
public:

    typedef ITERATOR const_iterator;

    CommPatternBase (MPI_Comm C)
: CommPattern (C)
    {
    };

    /// Return an iterator for the given rank (default=our rank)
    const_iterator begin (int rank = -1) const
    {
return const_iterator (*this, rank, false);
    }

    const_iterator end (int rank = -1) const
    {
return const_iterator (*this, rank, true);
    }
};

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif
