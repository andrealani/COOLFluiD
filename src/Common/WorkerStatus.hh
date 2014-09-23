// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_StatusWorker_hh
#define COOLFluiD_Common_StatusWorker_hh

//////////////////////////////////////////////////////////////////////////////

#include "EnumT.hh"
#include "Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
 
namespace Common {
 
//////////////////////////////////////////////////////////////////////////////

 class Common_API WorkerStatus
 {
  public:

  /// Enumeration of the worker statuses recognized in COOLFluiD
   enum Type  { INVALID = -1,
    STARTING    = 0,
    CONFIGURING = 1,
    EXITING     = 2,
    RUNNING     = 3,
    NOT_RUNNING = 4,
    WAITING     = 5,
    IDLE        = 6 };

    typedef Common::EnumT< WorkerStatus > ConverterBase;

    struct Common_API Convert : public ConverterBase
    {
    /// storage of the enum forward map
     static ConverterBase::FwdMap_t all_fwd;
    /// storage of the enum reverse map
     static ConverterBase::BwdMap_t all_rev;
    };

 }; // class StatusWorker

 ///////////////////////////////////////////////////////////////////////////

 Common_API std::ostream& operator<< ( std::ostream& os, const WorkerStatus::Type& in );
 Common_API std::istream& operator>> ( std::istream& is, WorkerStatus::Type& in );

 //////////////////////////////////////////////////////////////////////////////
} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_StatusWorker_hh
