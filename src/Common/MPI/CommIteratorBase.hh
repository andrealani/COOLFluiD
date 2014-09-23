// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_COMMITERATORBASE_HH
#define PARALLEL_MPI_COMMITERATORBASE_HH

#include "Common/MPI/CommPattern.hh"

namespace COOLFluiD  {
    namespace Common  {

//////////////////////////////////////////////////////////////////////////////

class Common_API CommIteratorSpecBase
{
public:
    typedef CommPattern::value_type value_type;

public:
    /// Return the rank in the communicator OF THE CommPattern for whom we
    /// are iterating the pattern.
    int GetMyRank () const {  return MyRank_; }

    /// Return the rank used in the actual communication.
    /// (This can be different from the rank we have in the CommPattern
    /// communicator)
    /// By default, we do not change the communicator so this returns the
    /// same as GetMyRank ()
    int GetCommunicationRank () const  { return MyRank_;  }

    int GetLogicalTime () const  {  return LogicalTime_;  }

    bool operator == (const CommIteratorSpecBase & other) const
    {
//      CFout << "CommIteratorSpecBase operator == " << "\n";
        cf_assert (&Pattern_ == &other.Pattern_);
        cf_assert (MyRank_ == other.MyRank_);

        // If they have different end state, they are not  alike
        if (LogicalTime_ * other.LogicalTime_ < 0)
            return false;

        // If we got here, they are either both ended or none is.
        if (LogicalTime_ < 0)
            return true;


        return (Current_ == other.Current_)
            && (LogicalTime_ == other.LogicalTime_);
    }

    /// By default, we do not change the communicator
    /// If we do, this function needs to be overrided
    /// This is the communicator that needs to be used for the actual
    /// communication.
    MPI_Comm GetCommunicator () const {  return Pattern_.GetCommunicator (); }


protected:
    CommIteratorSpecBase (const CommPattern & P, int MyRank)
        : Pattern_(P), MyRank_(MyRank), LogicalTime_(-1), Current_(-1,-1)
    {
        if (MyRank_ < 0)
            MyRank_ = Pattern_.GetRank ();
    }


    void CalcNext (bool First)
    {
        cf_assert (false);
    }

    const value_type & GetCurrent () const
    {
        return Current_;
    }

protected:
    const CommPattern & Pattern_;
    int MyRank_;
    int LogicalTime_;
    value_type Current_;
};

/// Helper class for communication iterators
///   Stores logical time, current value
///   Clients need to provide a base class that has CalcNext ()
///   (by deriving from SpecBase)
template <typename SPEC>
class CommIteratorBase :
    public SPEC,
    public std::iterator<std::input_iterator_tag, CommPattern::value_type>
{
public:

    typedef CommPattern::value_type value_type;

    CommIteratorBase (const CommPattern & P, int MyRank, bool End)
        : SPEC (P, MyRank)
    {
        if (!End)
            Step (true);
    }

    const typename SPEC::value_type operator * () const
    {
        return GetCurrent();
    }

    CommIteratorBase & operator ++ ()
    {
         Step (false);
         return *this;
    }

    CommIteratorBase & operator ++ (int dummy)
    {
        Step (false);
        return *this;
    }


    const typename SPEC::value_type * operator -> () const
    {
        return &GetCurrent();
    }

    typename SPEC::value_type * operator -> ()
    {
        return &GetCurrent();
    }

    bool operator == (const CommIteratorBase & other) const
    {
        return dynamic_cast<const SPEC &>(*this) ==
            dynamic_cast<const SPEC &>(other);

    }

    bool operator != (const CommIteratorBase & other) const
    {
        return ! (*this == other);
    }

protected:

    void Step (bool First)
    {
        SPEC::CalcNext (First);
    }

};

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif
