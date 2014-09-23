// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/MPI/MPIDataTypeRegistrar_Helper.hh"
#include "Common/MPI/MPIDataTypeHandler.hh"

namespace COOLFluiD {
   namespace Common  {

//////////////////////////////////////////////////////////////////////////////

MPIDataTypeRegistrar_Helper::MPIDataTypeRegistrar_Helper ()
	    : TheType (MPI_DATATYPE_NULL)
{
}

MPIDataTypeRegistrar_Helper::~MPIDataTypeRegistrar_Helper ()
{
}

void MPIDataTypeRegistrar_Helper::DoRegister ()
{
    MPIDataTypeHandler::GetHandler ().RegisterType (this);
}
//////////////////////////////////////////////////////////////////////////////
    }
}
