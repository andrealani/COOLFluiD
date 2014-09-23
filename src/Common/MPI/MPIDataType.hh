// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MPIDATATYPE_HH
#define MPIDATATYPE_HH

#include <mpi.h>

namespace COOLFluiD {
  namespace Common  {

//////////////////////////////////////////////////////////////////////////////

/// Interface for MPI Datatypes
class MPIDataType {
public:

  /// Regist the Datatype
  virtual void Register(MPI_Comm Comm) = 0;

  /// Unregist the Datatype
  virtual void UnRegister() = 0;

  /// Virtual destructor
  virtual ~MPIDataType ();

}; // end class MPIDataType

//////////////////////////////////////////////////////////////////////////////

  } // Common
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // MPIDATATYPE_HH
