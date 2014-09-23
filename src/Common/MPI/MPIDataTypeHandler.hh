// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MPIDATATYPEHANDLER_HH
#define MPIDATATYPEHANDLER_HH

/*****************************************************************
 *  This class handles correct initialisation and destruction    *
 *  of the registered datatypes                                  *
 *****************************************************************/

#include <mpi.h>
#include "Common/MPI/MPIDataTypeRegistrar.hh"
#include "Common/MPI/MPIDataType.hh"

#include <set>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/**
 *  This class handles correct initialisation and destruction
 *  of the registered datatypes
 */
class MPIDataTypeHandler {
public:

    static MPIDataTypeHandler & GetHandler ();

    template <class CTYPE>
    static MPI_Datatype GetType ();

    MPIDataTypeHandler (MPI_Comm Comm);
    ~MPIDataTypeHandler ();

    void RegisterType (MPIDataType * Type);

    bool IsInitialized () const;

    void InitTypes ();
    void DoneTypes ();

private:


    void InitType (MPIDataType * T);
    void DoneType (MPIDataType * T);


    MPI_Comm Communicator;
    bool Initialized;

    static MPIDataTypeHandler * TheInstance;

    typedef std::set<MPIDataType *> ContainerType;

    ContainerType Types;

}; // end class MPIDataTypeHandler

//////////////////////////////////////////////////////////////////////////////

template <class T>
MPI_Datatype MPIDataTypeHandler::GetType ()
{
    static MPIDataTypeRegistrar<T> TheRegistrar;
    return TheRegistrar.GetType ();
}

//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

#endif
