// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MPIDATATYPEREGISTRAR_HH
#define MPIDATATYPEREGISTRAR_HH

#include <mpi.h>

#include <typeinfo>

#include "Common/CFLog.hh"

/*
 * UNSAFE_MPITYPES enables generic MPI Type registration code.
 * This means MPI no longer 'knows' the exact meaning of the data,
 * but will consider it as meaningless bytes.
 *
 * This shouldn't be a problem when staying within one architecture within
 * the cluster...
 *
 */

#define UNSAFE_MPITYPES

#ifdef UNSAFE_MPITYPES
#include "Common/MPI/MPIDataTypeRegistrar_Helper.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
   namespace Common  {

//////////////////////////////////////////////////////////////////////////////

#ifndef UNSAFE_MPITYPES
//=========================== SAFE MPI Types ==========================
/*template <class T>
class MPIDataTypeRegistrar
{
public:

    // Every type should override this
    MPI_Datatype GetType () const
        { cf_assert (false); };

};*/

//===================================================================
#else
//=========================== UNSAFE MPI Types ======================

template <class T>
class MPIDataTypeRegistrar : public MPIDataTypeRegistrar_Helper {
public:

    MPIDataTypeRegistrar()
    {
      DoRegister ();
    }

    virtual ~MPIDataTypeRegistrar ()
    {
    }

    virtual void Register (MPI_Comm Comm)
    {
      CFLogDebugMin( "MPIDataTypeRegistrar(generic): Registering type " <<
        typeid(T).name() << " of size " << sizeof (T) << "\n");
      MPI_Type_contiguous (sizeof (T), MPI_BYTE, &TheType);
      MPI_Type_commit (&TheType);
    }

    virtual void UnRegister ()
    {
      CFLogDebugMin( "MPIDataTypeRegistrar(generic): "
        << "DeRegistering type " <<
        typeid(T).name() << " of size " << sizeof (T) << "\n");
      MPI_Type_free (&TheType);
    }

}; // end class MPIDataTypeRegistrar

#endif // UNSAFE_MPITYPES

#define REGISTRAR_TYPE(CT,MT) \
    template <>           \
    class MPIDataTypeRegistrar<CT>      \
    {             \
  public:           \
      MPI_Datatype GetType () const   \
    { return MT; };       \
    };


REGISTRAR_TYPE(double, MPI_DOUBLE)
REGISTRAR_TYPE(float,MPI_FLOAT)
REGISTRAR_TYPE(int,MPI_INT)
REGISTRAR_TYPE(unsigned int,MPI_UNSIGNED)
REGISTRAR_TYPE(long unsigned int,MPI_UNSIGNED_LONG)
REGISTRAR_TYPE(long long int, MPI_LONG_LONG_INT)

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

#endif
