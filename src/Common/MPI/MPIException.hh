// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_MPIEXCEPTION_HH
#define PARALLEL_MPI_MPIEXCEPTION_HH

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
   namespace Common  {

//////////////////////////////////////////////////////////////////////////////

class MPIException : public Common::Exception {
public:

  /// Constructor
  MPIException (const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "MPIException") {}

  /// Copy constructor
  MPIException(const MPIException& e) throw  () : Exception(e) {}

}; // class MPIException

//////////////////////////////////////////////////////////////////////////////

    } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // PARALLEL_MPI_MPIEXCEPTION_HH
