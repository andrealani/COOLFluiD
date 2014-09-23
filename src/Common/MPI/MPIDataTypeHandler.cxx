// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/MPI/MPIDataTypeHandler.hh"
#include "Common/CFLog.hh"
#include <mpi.h>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

MPIDataTypeHandler * MPIDataTypeHandler::TheInstance = 0;

void MPIDataTypeHandler::InitType (MPIDataType * Type)
{
    Type->Register (Communicator);
}

//////////////////////////////////////////////////////////////////////////////

void MPIDataTypeHandler::DoneType (MPIDataType * Type)
{
    Type->UnRegister ();
}

//////////////////////////////////////////////////////////////////////////////

void MPIDataTypeHandler::RegisterType (MPIDataType * Type)
{
    Types.insert(Type);
    if (IsInitialized ())
        InitType (Type);
}

void MPIDataTypeHandler::InitTypes ()
{
    cf_assert (!Initialized);
    ContainerType::const_iterator Iter = Types.begin ();
    while (Iter != Types.end ())
    {
        InitType (*Iter);
        Iter++;
    }
    Initialized = true;
}

//////////////////////////////////////////////////////////////////////////////

void MPIDataTypeHandler::DoneTypes()
{
  CFLogDebugMin("MPIDataTypeHandler::DoneTypes() begin" << "\n");

  cf_assert(Initialized);

  ContainerType::const_iterator Iter;
  for (Iter = Types.begin ();  Iter != Types.end (); Iter++) {
    DoneType(*Iter);
  }

  Types.clear();

  Initialized = false;

  CFLogDebugMin("MPIDataTypeHandler::DoneTypes() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

bool MPIDataTypeHandler::IsInitialized () const
{
  return Initialized;
}

//////////////////////////////////////////////////////////////////////////////

MPIDataTypeHandler::MPIDataTypeHandler (MPI_Comm Comm) :
    Communicator(Comm), Initialized(false)
{
  cf_assert (TheInstance == 0);
  TheInstance = this;
}

//////////////////////////////////////////////////////////////////////////////

MPIDataTypeHandler::~MPIDataTypeHandler()
{
  CFLogDebugMin( "MPIDataTypeHandler::~MPIDataTypeHandler (): " << Types.size () << " types to destroy" << "\n");

  if (IsInitialized ()) {
    DoneTypes ();
  }

  TheInstance = CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

MPIDataTypeHandler & MPIDataTypeHandler::GetHandler()
{
  cf_assert (TheInstance != CFNULL);
  return *TheInstance;
}

//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD
