// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_PATTERNRING_HH
#define PARALLEL_MPI_PATTERNRING_HH

#include <mpi.h>
#include "CommPattern.hh"
#include "CommIteratorBase.hh"

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

namespace
{

class PatternRingSpec : public CommIteratorSpecBase
{
protected:
    PatternRingSpec (const CommPattern & P, int MyRank)
	: CommIteratorSpecBase (P, MyRank)
    {
	// TODO: create ring communicator

	// Receive from
	Current_.first=MyRank_-1;
	// send to
	Current_.second=MyRank_+1;

	if (Current_.first < 0)
	    Current_.first = Pattern_.GetCommSize()-1;

	if (Current_.second >= Pattern_.GetCommSize())
	    Current_.second = 0;
    }

    void CalcNext (bool first)
    {
	cf_assert (first || LogicalTime_ >= 0);

	if (first)
	    LogicalTime_ = 0;

	++LogicalTime_;
	if (LogicalTime_ == Pattern_.GetCommSize())
	    LogicalTime_ = -1;
    }
};

typedef CommIteratorBase<PatternRingSpec> PatternRingIterator;
}

typedef CommPatternBase<PatternRingIterator> PatternRing;


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
	}
    }
}

#endif
