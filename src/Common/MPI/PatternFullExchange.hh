// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_PATTERNFULLEXCHANGE_HH
#define PARALLEL_MPI_PATTERNFULLEXCHANGE_HH

#include "CommIteratorBase.hh"
#include "CommPattern.hh"

using namespace std;

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////


/**
 * Full exchange pattern (assumed to be bidirectional, to create an
 * unidirectional iterator an adapter could be used)
 *
 * Algorithm found by Jesse & Tom
 */
class FullExchangeSpec : public CommIteratorSpecBase
{

protected:
    FullExchangeSpec (const CommPattern & P, int MyRank)
        : CommIteratorSpecBase (P, MyRank),
          CommSize_(Pattern_.GetCommSize()),
          CustomCommSize_(CommSize_+(CommSize_ & 1)),
          RoundCount_(CustomCommSize_-1)
    {
        LogicalTime_ = -1;
    }

    void CalcNext (bool First)
    {
        //CFout << "First=" << First << ", logicaltime=" << LogicalTime_ << "\n";
        cf_assert (First || LogicalTime_ >= 0);

        do
        {
            ++LogicalTime_;

            if ((LogicalTime_) >= RoundCount_)
            {
                LogicalTime_ = -1;
                return ;
            }

            const int t = LogicalTime_ + 1;
            int p = -1;
            if (MyRank_==(CustomCommSize_-1))
            {
                if (((CustomCommSize_-1)+t)%2)
                {
                    // oneven
                    p = (((CustomCommSize_-1)+t) % (CustomCommSize_-1))/2;
                }
                else
                {
                    p = (((CustomCommSize_-1)+t)/2) % (CustomCommSize_-1);
                }
            }
            else
            {
                p = ((CustomCommSize_-1)-MyRank_+t) % (CustomCommSize_-1);
                if (p==MyRank_)
                    p=CustomCommSize_-1;

            }

            if (p >= CommSize_)
                p = -1;
/*          // Now order so we always send to the higher rank
            if (p < MyRank_)
                Current_ = std::make_pair(p, -1);
            else
                Current_ = std::make_pair(-1, p); */

            Current_ = std::make_pair(p, p);

        } while (false); /*while (Current_.first >= CommSize_ || Current_.second >= CommSize_);*/
    }

private:
    const int CommSize_;
    const int CustomCommSize_;
    const int RoundCount_;
};


typedef CommIteratorBase<FullExchangeSpec> FullExchangeIterator;
typedef CommPatternBase<FullExchangeIterator> PatternFullExchange;

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif
