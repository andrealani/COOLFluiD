// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_MPI_MPIStructDef_hh
#define COOLFluiD_Common_MPI_MPIStructDef_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// The following struct defines the basic data needed by a generic
/// MPI struct
/// @author Andrea Lani
struct MPIStruct {
  MPI_Datatype type;
  void* start;
};

/// The following class defines static functions that build
/// derived MPI types and set MPIStruct objects
/// @author Andrea Lani
class MPIStructDef {
public:
  
  /// The following function returns the MPI basic type
  /// corresponding to the given argument
#define MPIDTYPE(__type__,__mpidtype__)  	\
static MPI_Datatype getMPIType(__type__* type)  \
{  				  \
  return __mpidtype__;  		  \
}

MPIDTYPE(double,MPI_DOUBLE)
MPIDTYPE(float,MPI_FLOAT)
MPIDTYPE(int,MPI_INT)
MPIDTYPE(long,MPI_LONG)
MPIDTYPE(short,MPI_SHORT)
MPIDTYPE(unsigned int,MPI_UNSIGNED)
MPIDTYPE(long unsigned int,MPI_UNSIGNED_LONG)
MPIDTYPE(short unsigned int,MPI_UNSIGNED_SHORT)
MPIDTYPE(char,MPI_CHAR)

#undef MPIDTYPE

  /// The following function build a MPIStruct with 2 types
  template <typename T1, typename T2>
  static void buildMPIStruct(T1* t1, T2* t2,
  		     int blockLengths[],
  		     MPIStruct& obj)
  {
    const unsigned int N = 2;
    MPI_Datatype typelist[N];
    typelist[0] = getMPIType(t1);
    typelist[1] = getMPIType(t2);

    MPI_Aint displacements[N];
    MPI_Aint startAddress;
    MPI_Aint address;
    displacements[0] = 0;

    MPI_Address(t1,&startAddress);
    MPI_Address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Type_struct(N, blockLengths, displacements, typelist, &obj.type);
    MPI_Type_commit(&obj.type);
    obj.start = t1;
  }

  /// The following function build a MPIStruct with 3 types
  template <typename T1, typename T2, typename T3>
  static void buildMPIStruct(T1* t1, T2* t2, T3* t3,
  		     int blockLengths[],
  		     MPIStruct& obj)
  {
    const unsigned int N = 3;
    MPI_Datatype typelist[N];
    typelist[0] = getMPIType(t1);
    typelist[1] = getMPIType(t2);
    typelist[2] = getMPIType(t3);

    MPI_Aint displacements[N];
    MPI_Aint startAddress;
    MPI_Aint address;
    displacements[0] = 0;

    MPI_Address(t1,&startAddress);
    MPI_Address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Type_struct(N, blockLengths, displacements, typelist, &obj.type);
    MPI_Type_commit(&obj.type);
    obj.start = t1;
  }

  /// The following function build a MPIStruct with 4 types
  template <typename T1, typename T2, typename T3, typename T4>
  static void buildMPIStruct(T1* t1, T2* t2, T3* t3, T3* t4,
  		     int blockLengths[],
  		     MPIStruct& obj)
  {
    const unsigned int N = 4;
    MPI_Datatype typelist[N];
    typelist[0] = getMPIType(t1);
    typelist[1] = getMPIType(t2);
    typelist[2] = getMPIType(t3);
    typelist[3] = getMPIType(t4);

    MPI_Aint displacements[N];
    MPI_Aint startAddress;
    MPI_Aint address;
    displacements[0] = 0;

    MPI_Address(t1,&startAddress);
    MPI_Address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Type_struct(N, blockLengths, displacements, typelist, &obj.type);
    MPI_Type_commit(&obj.type);
    obj.start = t1;
  }

  /// The following function build a MPIStruct with 5 types
  template <typename T1, typename T2,
    typename T3, typename T4,
    typename T5>
  static void buildMPIStruct(T1* t1, T2* t2, T3* t3,
  		     T4* t4, T5* t5,
  		     int blockLengths[],
  		     MPIStruct& obj)
  {
    const unsigned int N = 5;
    MPI_Datatype typelist[N];
    typelist[0] = getMPIType(t1);
    typelist[1] = getMPIType(t2);
    typelist[2] = getMPIType(t3);
    typelist[3] = getMPIType(t4);
    typelist[4] = getMPIType(t5);

    MPI_Aint displacements[N];
    MPI_Aint startAddress;
    MPI_Aint address;
    displacements[0] = 0;

    MPI_Address(t1,&startAddress);
    MPI_Address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Address(t5,&address);
    displacements[4] = address - startAddress;

    MPI_Type_struct(N, blockLengths, displacements, typelist, &obj.type);
    MPI_Type_commit(&obj.type);
    obj.start = t1;
  }

  /// The following function build a MPIStruct with 6 types
  template <typename T1, typename T2,
    typename T3, typename T4,
    typename T5, typename T6>
  static void buildMPIStruct(T1* t1, T2* t2, T3* t3,
  		     T4* t4, T5* t5, T6* t6,
  		     int blockLengths[],
  		     MPIStruct& obj)
  {
    const unsigned int N = 6;
    MPI_Datatype typelist[N];
    typelist[0] = getMPIType(t1);
    typelist[1] = getMPIType(t2);
    typelist[2] = getMPIType(t3);
    typelist[3] = getMPIType(t4);
    typelist[4] = getMPIType(t5);
    typelist[5] = getMPIType(t6);

    MPI_Aint displacements[N];
    MPI_Aint startAddress;
    MPI_Aint address;
    displacements[0] = 0;

    MPI_Address(t1,&startAddress);
    MPI_Address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Address(t5,&address);
    displacements[4] = address - startAddress;

    MPI_Address(t6,&address);
    displacements[5] = address - startAddress;

    MPI_Type_struct(N, blockLengths, displacements, typelist, &obj.type);
    MPI_Type_commit(&obj.type);
    obj.start = t1;
  }

private:

  /// Private constructor
  MPIStructDef();

};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_MPI_MPIStructDef_hh
