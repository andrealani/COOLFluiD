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
  
  /// Get the MPI_Datatype corresponding to MPI_Offset
  static MPI_Datatype getMPIOffsetType() 
  {
    return MPI_OFFSET;
    // MPI_Datatype offsetType = MPI_DATATYPE_NULL;
    // if      (sizeof(MPI_Offset) == sizeof(int))  { return MPI_INT; }
    // else if (sizeof(MPI_Offset) == sizeof(long)) { return MPI_LONG; }
    // else if (sizeof(MPI_Offset) == sizeof(long long)) { return MPI_LONG_LONG; }
    // else { MPI_Abort(MPI_COMM_WORLD, 1); }
    // return  MPI_DATATYPE_NULL;
  }
  
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
MPIDTYPE(unsigned int,MPI_UNSIGNED)
MPIDTYPE(long unsigned int,MPI_UNSIGNED_LONG)
MPIDTYPE(long int,MPI_LONG)
MPIDTYPE(long long int, MPI_LONG_LONG_INT)
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

    MPI_Get_address(t1,&startAddress);
    MPI_Get_address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Type_create_struct(N, blockLengths, displacements, typelist, &obj.type);
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

    MPI_Get_address(t1,&startAddress);
    MPI_Get_address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Get_address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Type_create_struct(N, blockLengths, displacements, typelist, &obj.type);
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

    MPI_Get_address(t1,&startAddress);
    MPI_Get_address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Get_address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Get_address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Type_create_struct(N, blockLengths, displacements, typelist, &obj.type);
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

    MPI_Get_address(t1,&startAddress);
    MPI_Get_address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Get_address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Get_address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Get_address(t5,&address);
    displacements[4] = address - startAddress;

    MPI_Type_create_struct(N, blockLengths, displacements, typelist, &obj.type);
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

    MPI_Get_address(t1,&startAddress);
    MPI_Get_address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Get_address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Get_address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Get_address(t5,&address);
    displacements[4] = address - startAddress;

    MPI_Get_address(t6,&address);
    displacements[5] = address - startAddress;
    
    MPI_Type_create_struct(N, blockLengths, displacements, typelist, &obj.type);
    MPI_Type_commit(&obj.type);
    obj.start = t1;
  }

  /// The following function build a MPIStruct with 14 types
  template <typename T1, typename T2,
    typename T3, typename T4,
    typename T5, typename T6, typename T7, typename T8, typename T9,typename T10,
            typename T11, typename T12, typename T13, typename T14>
  static void buildMPIStruct(T1* t1, T2* t2, T3* t3,
             T4* t4, T5* t5, T6* t6, T7* t7, T8* t8,
             T9* t9, T10* t10, T11* t11, T12* t12,
             T13* t13, T14* t14,

             int blockLengths[],
             MPIStruct& obj)
  {
    const unsigned int N = 14;
    MPI_Datatype typelist[N];
    typelist[0] = getMPIType(t1);
    typelist[1] = getMPIType(t2);
    typelist[2] = getMPIType(t3);
    typelist[3] = getMPIType(t4);
    typelist[4] = getMPIType(t5);
    typelist[5] = getMPIType(t6);
    typelist[6] = getMPIType(t7);
    typelist[7] = getMPIType(t8);
    typelist[8] = getMPIType(t9);
    typelist[9] = getMPIType(t10);
    typelist[10] = getMPIType(t11);
    typelist[11] = getMPIType(t12);
    typelist[12] = getMPIType(t13);
    typelist[13] = getMPIType(t14);

    MPI_Aint displacements[N];
    MPI_Aint startAddress;
    MPI_Aint address;
    displacements[0] = 0;

    MPI_Get_address(t1,&startAddress);
    MPI_Get_address(t2,&address);
    displacements[1] = address - startAddress;

    MPI_Get_address(t3,&address);
    displacements[2] = address - startAddress;

    MPI_Get_address(t4,&address);
    displacements[3] = address - startAddress;

    MPI_Get_address(t5,&address);
    displacements[4] = address - startAddress;

    MPI_Get_address(t6,&address);
    displacements[5] = address - startAddress;

    MPI_Get_address(t7,&address);
    displacements[6] = address - startAddress;

    MPI_Get_address(t8,&address);
    displacements[7] = address - startAddress;

    MPI_Get_address(t9,&address);
    displacements[8] = address - startAddress;

    MPI_Get_address(t10,&address);
    displacements[9] = address - startAddress;

    MPI_Get_address(t11,&address);
    displacements[10] = address - startAddress;

    MPI_Get_address(t12,&address);
    displacements[11] = address - startAddress;

    MPI_Get_address(t13,&address);
    displacements[12] = address - startAddress;

    MPI_Get_address(t14,&address);
    displacements[13] = address - startAddress;

    MPI_Type_create_struct(N, blockLengths, displacements, typelist, &obj.type);
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
