// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalReduce_hh
#define COOLFluiD_Framework_GlobalReduce_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/PM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Empty template class, to be specialized for the correct PE
template <typename RESULTTYPE, typename PEMODE = Common::PM_CUR>
class GlobalReduce
{
};

//////////////////////////////////////////////////////////////////////////////

/// Tag class for GlobalReduceOperation
class Common_API GRO_SUM {};

/// Tag class for GlobalReduceOperation
class Common_API GRO_MAX {};

/// Tag class for GlobalReduceOperation
class Common_API GRO_MIN {};

/// Helper class
template <typename TAGCLASS, typename BASETYPE,typename PEMODE = Common::PM_CUR>
class GlobalReduceOperationHelper;

/// TODO: helper function to avoid specifying BASETYPE
template <typename TAGCLASS, typename BASETYPE>
inline void GlobalReduceOperation (BASETYPE * Source, BASETYPE * Dest,
unsigned int Count = 1)
{
  GlobalReduceOperationHelper<TAGCLASS, BASETYPE>::DoReduce(Source, Dest, Count);
}

//////////////////////////////////////////////////////////////////////////////
    
  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
#  include "Framework/GlobalReduceMPI.hh"
#else
#  include "Framework/GlobalReduceSERIAL.hh"
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GlobalReduce_hh
