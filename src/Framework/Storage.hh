// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Storage_hh
#define COOLFluiD_Framework_Storage_hh

//////////////////////////////////////////////////////////////////////////////

/// This is the header to be included when dealing with DataSorage
/// and DataHandle', not the class headers. This header takes care of
/// selecting the appropriate implementation, with LOCAL or GLOBAL.
/// If the communication type is not defined, LOCAL is assumed.

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalCommTypes.hh"

// The local storage type is always present
// and is given as default
#include "Framework/DataStorage.hh"
#include "Framework/DataHandle.hh"

#ifdef CF_HAVE_MPI
# include "Framework/DataHandleMPI.hh"
#endif // CF_HAVE_MPI

#ifdef CF_HAVE_SHM
/// @todo support for shared memory is not implemented
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Storage_hh

