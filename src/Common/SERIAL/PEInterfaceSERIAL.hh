// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PEINTERFACESERIAL_HH
#define COOLFluiD_Common_PEINTERFACESERIAL_HH

#include "Common/PEInterface.hh"
#include "Common/NonCopyable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// The PE for the 'SERIAL' parallel model.
/// We inherit from PEInterface<bool> just to have
/// to be forced to implement the base class functions
/// It cannot be in the SERIAL namespace because it is a specialisation of a
/// class of the Parallel namespace
template <>
class PEInterface<PM_SERIAL> : public PEInterfaceBase,
                             public Common::NonCopyable<PEInterface<PM_SERIAL> >
{
public:

    PEInterface (int * argc, char *** args);

    ~PEInterface ();

    inline unsigned int GetProcessorCount () const;
    inline unsigned int GetRank () const;
    inline bool IsParallel () const;
    inline std::string GetName () const;
};

//////////////////////////////////////////////////////////////////////////////

std::string PEInterface<PM_SERIAL>::GetName () const
{
    return "SERIAL";
}

unsigned int PEInterface<PM_SERIAL>::GetProcessorCount () const
{
    return 1;
}

unsigned int PEInterface<PM_SERIAL>::GetRank () const
{
    return 0;
}

bool PEInterface<PM_SERIAL>::IsParallel () const
{
    return false;
}

//////////////////////////////////////////////////////////////////////////////

    } // Common

} // COOLFluiD

#endif // COOLFluiD_Common_PEINTERFACESERIAL_HH
