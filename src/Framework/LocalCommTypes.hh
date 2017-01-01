// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LocalCommTypes_HH
#define COOLFluiD_Framework_LocalCommTypes_HH

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/GrowArray.hh"
#include "Common/IsFundamental.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaVector.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is the policy for local simulations with no
/// inter-process communication and storage based on STL vectors
class LocalStlVector
{
  public:
  template < typename ETYPE >
  class StoragePolicy
  {
    public:
      typedef std::vector<ETYPE> ContainerType;
      typedef ETYPE              ElemType;
  };
};

//////////////////////////////////////////////////////////////////////////////

/// This class is the policy for local simulations with no
/// inter-process communication and storage based on GrowArray
class LocalGrowArray
{
  public:
  template < typename ETYPE >
  class StoragePolicy
  {
    public:
      typedef Common::GrowArray<ETYPE> ContainerType;
      typedef ETYPE                    ElemType;
  };
};

//////////////////////////////////////////////////////////////////////////////x

#ifdef CF_HAVE_CUDA

template <typename T, bool isArithmetic = false> 
struct SelectCudaArray {
#ifdef CF_ENABLE_GROWARRAY
  typedef  Common::GrowArray<T> ARRAY;
#else
  typedef  std::vector<T> ARRAY;
#endif
};

template <typename T> 
struct SelectCudaArray<T, true> {typedef CudaEnv::CudaVector<T> ARRAY;};

/// This class is the policy for local simulations with no
/// inter-process communication and storage based on GrowArray
class LocalCudaArray
{
public:
  template < typename ETYPE >
  class StoragePolicy
  {
  public:
    typedef typename SelectCudaArray<ETYPE, Common::IsArithmetic<ETYPE>::VALUE>::ARRAY ContainerType;
    typedef ETYPE ElemType;
  };
};
#endif

//////////////////////////////////////////////////////////////////////////////

/// If GrowArray is supported set it to the default LOCAL communicator type
/// else is based on std:vector
#ifdef CF_HAVE_CUDA
 typedef LocalCudaArray LOCAL;
#else
#ifdef CF_ENABLE_GROWARRAY
typedef LocalGrowArray LOCAL;
#else
typedef LocalStlVector LOCAL;
#endif
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA
#define HPCVEC(stype) CudaEnv::CudaVector<stype>
#define HPCVECMALLOC(stype) CudaEnv::CudaVector<stype, CudaEnv::MallocHostAlloc<stype> >
#else
#define HPCVEC(stype) std::vector<stype>
#define HPCVECMALLOC(stype) std::vector<stype>
#endif

//////////////////////////////////////////////////////////////////////////////
    
template <typename ETYPE>
class LocalArray {
public:
#ifdef CF_HAVE_CUDA
#ifndef CF_HAVE_CUDA_MALLOC
  typedef CudaEnv::CudaVector<ETYPE> TYPE;
#else
  typedef CudaEnv::CudaVector<ETYPE, CudaEnv::MallocHostAlloc<ETYPE> > TYPE;
#endif
  typedef CudaEnv::CudaVector<ETYPE, CudaEnv::MallocHostAlloc<ETYPE> > MALLOC_TYPE;
#else
  typedef std::vector<ETYPE> TYPE;
  typedef std::vector<ETYPE> MALLOC_TYPE;
#endif
};

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} //  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LocalCommTypes_HH
