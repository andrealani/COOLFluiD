// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalCommTypes_HH
#define COOLFluiD_Framework_GlobalCommTypes_HH

//////////////////////////////////////////////////////////////////////////////

#include "Framework/LocalCommTypes.hh"

#ifdef CF_HAVE_MPI
# include "Common/MPI/ParVector.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
  /// This class is the policy for local simulations with MPI based
  /// inter-process communication and storage based on MPI vectors
  class GlobalMpiVector
  {
    public:
    template < typename ETYPE >
    class StoragePolicy
    {
    public:
      typedef COOLFluiD::Common::ParVector<ETYPE> ContainerType;
      typedef ETYPE ElemType;
    };
  };

  /// this is the default MPI interface
  typedef GlobalMpiVector MPI;

  /// With MPI support, GLOBAL is set to MPI
  typedef MPI GLOBAL;
#else
  /// Without MPI support, GLOBAL is set to LOCAL
  typedef LOCAL GLOBAL;
#define CF_GLOBAL_EQUAL_LOCAL
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} //  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GlobalCommTypes_HH
