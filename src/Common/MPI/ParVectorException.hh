// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PAR_VECTOR_MPI_EXCEPTION_HH
#define PAR_VECTOR_MPI_EXCEPTION_HH

#include "Common/ParallelException.hh"

namespace COOLFluiD  {
    namespace Common  {

class ParVectorException : public ParallelException
{
    public:
	ParVectorException (const Common::CodeLocation& where, const char * S) : ParallelException(where, S)
	{
	}
};

class DoubleElementException : public ParVectorException
{
public:
    DoubleElementException (const Common::CodeLocation& where, const char * S) : ParVectorException (where,S)
    {
    }
};

class StorageException : public ParVectorException
{
public:
    StorageException (const Common::CodeLocation& where, const char * S) : ParVectorException (where, S)
    {
    }
};

class NotFoundException : public ParVectorException
{

public:
    NotFoundException (const Common::CodeLocation& where, const char * S) : ParVectorException (where, S)
    {
    }
};

} // Common

} // COOLFluiD


#endif
