// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_CudaVector_hh 
#define COOLFluiD_CudaEnv_CudaVector_hh

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include "CudaEnv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
          
//////////////////////////////////////////////////////////////////////////////

/// This class implements a CUDA stream
///
/// @author Andrea Lani
struct CudaStream {
  
  // CUDA stream object
  cudaStream_t stream;
  
  // size of the corresponding slice array of data
  size_t       size;
  
  // pointer to the start of the corresponding slice array
  size_t       start;
};

//////////////////////////////////////////////////////////////////////////////

enum TransferType {TO_GPU = 0, FROM_GPU = 1};
template <typename ARRAY, TransferType TT> class CudaCopy;

//////////////////////////////////////////////////////////////////////////////

/// This class implements a basic allocator for CudaVector.
///
/// @author Andrea Lani
template <typename T>
class CudaBaseAlloc {
public:
  
  /// constructor
  CudaBaseAlloc() : m_isAlloc(false), m_size(0) {}
  
  /// copy constructor
  CudaBaseAlloc(const CudaBaseAlloc& in) : m_isAlloc(in.m_isAlloc), m_size(in.m_size),  m_dataDev(NULL) {}
  
  /// overloading of assignment operator
  const CudaBaseAlloc& operator=(const CudaBaseAlloc& in) 
  {
    m_size = in.m_size;
    return *this;
  }
  
  /// destructor
  virtual ~CudaBaseAlloc() {}
  
  /// allocate the array on the device
  void allocateDev()
  { 
    freeDev();
    // data size allocated for the GPU corresponds to maximum slice size 
    CudaEnv::allocDev(m_dataDev, m_size);
    m_isAlloc = true;
  } 
  
  /// deallocate the array on the device
  void freeDev()
  { 
    if (m_isAlloc) {
      // data size allocated for the GPU corresponds to maximum slice size 
      CudaEnv::free(m_dataDev); 
      m_isAlloc = false;
    }
  }
  
  /// get pointer to the GPU array
  /// @post the device is allocated on-the-fly the first time you call this function
  T* ptrDev() {if (!m_isAlloc) {allocateDev();} cfassert(m_isAlloc); return m_dataDev;}
  
protected: // data
  
  /// flag telling if the array has been allocated
  bool m_isAlloc;
  
  /// size of the array on CPU
  size_t m_size;
  
  /// data array on the GPU
  T* m_dataDev;
};

//////////////////////////////////////////////////////////////////////////////

/// This class implements an allocator using pinned memory on the host for 
/// CudaVector.
///
/// @author Andrea Lani 
template <typename T>
class PinnedHostAlloc : public CudaBaseAlloc<T> {
public:
  
  /// constructor
  PinnedHostAlloc() : CudaBaseAlloc<T>(), m_streams(), m_data(NULL) {}
  
  /// copy constructor
  PinnedHostAlloc(const PinnedHostAlloc& in) : CudaBaseAlloc<T>(in), m_streams(in.m_streams)
  {
    CudaEnv::allocHost(m_data, this->m_size);
    this->allocateDev();
  }
  
  /// copy constructor
  const PinnedHostAlloc& operator=(const PinnedHostAlloc& in) 
  {
    freeHost();
    this->freeDev();
    this->m_size = in.m_size;
    m_streams.resize(in.m_streams.size());
    m_streams = in.m_streams;
    CudaEnv::allocHost(m_data, this->m_size);
    this->allocateDev();
    return *this;
  }
  
  /// destructor
  ~PinnedHostAlloc() {freeHost(); this->freeDev();}
  
  /// allocate 
  /// @param ns       full size on the CPU
  /// @param streams  stream array for asynchronous transfer/computation
  /// @param ssizes   sizes of the slices to be put/get on GPU
  void allocate(size_t ns,  
		std::vector<cudaStream_t>* streams, 
		std::vector<size_t>* ssizes)
  { 
    freeHost();
    
    // allocate data on the CPU (full size)  
    this->m_size = ns;
    CudaEnv::allocHost(m_data, this->m_size);
    
    // allocate data on the GPU with sizes depending on the stream data
    if (streams != NULL) {
      m_streams.resize(streams->size());
      for (size_t i = 0; i < streams->size(); ++i) {
	m_streams[i].stream  = (*streams)[i];
      }
      
      // if ssizes == NULL then a default size initialization is performed
      cfassert((streams != NULL && ssizes == NULL) || (streams->size() == ssizes->size()));
      for (size_t i = 0; i < streams->size(); ++i) {
	m_streams[i].size = (ssizes != NULL) ? (*ssizes)[i] : ns/streams->size();
      }
      if (ssizes != NULL) {m_streams.back().size += ns%streams->size();}
    }
    else {
      m_streams.resize(1);  m_streams[0].size = ns;
    }
    
    size_t sum = 0;
    for (size_t i = 0; i < m_streams.size(); ++i) {
      m_streams[i].start = sum;
      sum += m_streams[i].size;
    }
  }
    
  /// size of array on CPU
  size_t size() const {return this->m_size;}
  
  /// size of the slice corresponding to the given streamID 
  size_t size(size_t streamID) const 
  {
    cfassert(streamID < m_streams.size());
    return m_streams[streamID].size;
  }
    
  /// start of the slice corresponding to the given streamID 
  size_t start(size_t streamID) const 
  {
    cfassert(streamID < m_streams.size());
    return m_streams[streamID].start;
  }
  
  /// start of the slice corresponding to the given streamID 
  cudaStream_t* stream(size_t streamID)
  {
    // this has to be safe when a single stream is used
    cfassert(streamID < m_streams.size());
    return (m_streams.size() == 1) ? NULL : &(m_streams[streamID].stream);
  }
  
  /// subscripting operator for CPU array (returning value)
  T getv(size_t i) const {cfassert(i < this->m_size); return m_data[i];}
  
  /// subscripting operator for CPU array (returning reference)
  T& getr(size_t i) {cfassert(i < this->m_size); return m_data[i];}
  
  /// get pointer to the CPU array
  T* ptr() {return m_data;}
  
  /// get pointer to the CPU array entry corresponding to the given streamID slice
  T* ptr(size_t streamID) 
  {
    cfassert(start(streamID) < this->m_size); 
    return &m_data[start(streamID)];
  }
  
  /// get pointer to the CPU array entry corresponding to the given streamID slice
  /// @post the device is allocated on-the-fly the first time you call this function
  T* ptrDevS(size_t streamID)
  {
    if (!this->m_isAlloc) {this->allocateDev();}
    cfassert(start(streamID) < this->m_size); 
    return &this->m_dataDev[start(streamID)];
  }
  
  /// deallocate  the memory on CPU
  void freeHost() {if (this->m_size > 0) {CudaEnv::freeHost(m_data);} this->m_size = 0; m_data = NULL;}
  
protected:
  
  /// data array on the CPU
  T* m_data;
  
  /// sizes of the slices to be put/get on GPU 
  std::vector<CudaStream> m_streams;
};


//////////////////////////////////////////////////////////////////////////////

/// This class implements an allocator using malloc for the host in a CudaVector.
///
/// @author Andrea Lani 
template <typename T>
class MallocHostAlloc : public CudaBaseAlloc<T> {
public:
  
  /// constructor
  MallocHostAlloc() : CudaBaseAlloc<T>(), m_data() {}
  
  /// copy constructor
  MallocHostAlloc(const MallocHostAlloc& in) : CudaBaseAlloc<T>(in)
  {
    m_data.resize(this->m_size);
    this->allocateDev();
  }
  
  /// copy constructor
  const MallocHostAlloc& operator=(const MallocHostAlloc& in) 
  {
    freeHost();
    this->freeDev();
    this->m_size = in.m_size;
    m_data.resize(this->m_size);
    this->allocateDev();
    return *this;
  }
  
  /// destructor
  ~MallocHostAlloc() {freeHost(); this->freeDev();}
  
  /// allocate 
  /// @param ns       full size on the CPU
  /// @param streams  stream array for asynchronous transfer/computation
  /// @param ssizes   sizes of the slices to be put/get on GPU
  void allocate(size_t ns,  
		std::vector<cudaStream_t>* streams, 
		std::vector<size_t>* ssizes)
  { 
    freeHost();
    
    // allocate data on the CPU (full size)  
    this->m_size = ns;
    m_data.resize(this->m_size);
  }
  
  /// deallocate  the memory on both CPU and GPU 
  void freeHost() {if (this->m_size > 0) {std::vector<T>().swap(m_data);} this->m_size = 0;}
  
  /// size of array on CPU
  size_t size() const {return this->m_size;}
  
  /// size of the slice corresponding to the given streamID 
  size_t size(size_t streamID) const {return this->m_size;}
    
  /// start of the slice corresponding to the given streamID 
  size_t start(size_t streamID) const {return 0;}
  
  /// start of the slice corresponding to the given streamID 
  cudaStream_t* stream(size_t streamID) {return NULL;}
  
  /// subscripting operator for CPU array (returning value)
  T getv(size_t i) const {cfassert(i < this->m_size); return m_data[i];}
  
  /// subscripting operator for CPU array (returning reference)
  T& getr(size_t i) {cfassert(i < this->m_size); return m_data[i];}
  
  /// get pointer to the CPU array
  T* ptr() {return &m_data[0];}
  
  /// get pointer to the CPU array entry corresponding to the given streamID slice
  T* ptr(size_t streamID) {return ptr();}
  
  /// get pointer to the CPU array entry corresponding to the given streamID slice
  /// @post the device is allocated on-the-fly the first time you call this function
  T* ptrDevS(size_t streamID) {if (!this->m_isAlloc) {this->allocateDev();} return this->m_dataDev;}
  
protected:
  
  /// data array on the CPU
  std::vector<T> m_data;
  
};

//////////////////////////////////////////////////////////////////////////////

/// This class implements a slice array for a CudaVector.
///
/// @author Andrea Lani 
template <typename T, template <typename T1 = T> class ALLOC = PinnedHostAlloc> 
class CudaVectorSlice {
  
  typedef CudaVectorSlice<T, ALLOC> SELF;   
  
public:
  
  /// constructor
  CudaVectorSlice(size_t ns, T* ptr, T* ptrDev, cudaStream_t* cs) : 
    m_size(ns), m_ptr(ptr), m_ptrDev(ptrDev), m_cs(cs) {}
  
  /// destructor
  ~CudaVectorSlice() {}
  
  /// copy constructor
  CudaVectorSlice(const CudaVectorSlice& in) :
    m_size(in.m_size), m_ptr(in.m_ptr), m_ptrDev(in.m_ptrDev), m_cs(in.m_cs) {}
  
  /// (a)synchronous copy of all data from CPU to GPU
  void put() {CudaCopy<SELF, TO_GPU>(m_ptrDev, m_ptr, m_size, m_cs);}
  
  /// (a)synchronous copy of all data from GPU to CPU
  void get() {CudaCopy<SELF, FROM_GPU>(m_ptr, m_ptrDev, m_size, m_cs);}
  
  /// get the raw pointer on CPU
  T* ptr() {return m_ptr;}
  
  /// get the raw pointer on GPU
  T* ptrDev() {return m_ptrDev;}
  
  /// size
  size_t size() const {return m_size;}
  
private: /// helper functions
  
  /// default constructor is forbidden
  CudaVectorSlice();
  
  /// assignment operator is forbidden
  const CudaVectorSlice& operator=(const CudaVectorSlice& c);
  
private:
  
  /// size
  size_t m_size;
  
  /// pointer to array on CPU
  T* m_ptr;
  
  /// pointer to array on GPU
  T* m_ptrDev;  
  
  /// cuda stream
  cudaStream_t* m_cs;
  
};

//////////////////////////////////////////////////////////////////////////////

/// This class implements a CudaVector with default pinned allocator.
///
/// @author Andrea Lani  
template <typename T, template <typename T1 = T> class ALLOC = PinnedHostAlloc> 
class CudaVector {
  
  typedef CudaVector<T, ALLOC> SELF;   

public:
  typedef T TYPE;
  
  /// default constructor
  CudaVector() : m_alloc(), m_stride(1) {}
  
  /// Constructor
  /// @param init     value to use for initialization
  /// @param ns       full size on the CPU
  /// @param streams  stream array for asynchronous transfer/computation
  /// @param ssizes   sizes of the slices to be put/get on GPU
  CudaVector(T init,
	     size_t ns, 
	     std::vector<cudaStream_t>* streams = NULL, 
	     std::vector<size_t>* ssizes = NULL) : m_alloc(), m_stride(1)
  {
    resize(init, ns, streams, ssizes);
  }
  
  /// Dummy constructor for compatibility with @see ParVector
  /// @param init     value to use for initialization
  /// @param ns       full size on the CPU
  /// @param dummy    dummy integer
  CudaVector(T init, size_t ns, size_t stride) : 
    m_alloc(), m_stride(std::max(stride,(size_t)1))
  {
    cfassert(ns == 0);
  }
  
  /// copy constructor
  CudaVector (const CudaVector& in) : m_alloc(in.m_alloc), m_stride(in.m_stride) {}
  
  /// assignment operator (from input array)
  template <typename INPUT>
  const CudaVector& operator=(const INPUT& in) 
  {   
    m_alloc = in.m_alloc(); 
    return *this;
  }
  
  /// assignment operator (from input value)
  const CudaVector& operator=(const T& value) 
  {   
    for (CFuint i = 0; i < m_alloc.size(); ++i) {m_alloc.getr(i) = value;}
    return *this;
  }
  
  /// destructor
  ~CudaVector() {}
  
  /// this is just kept momentarily for compatibility with ParVector 
  void grow() {}
  
  /// free the memory
  void free() {if (size() > 0) {m_alloc.freeHost(); m_alloc.freeDev();}}
  
  /// resize storage (on CPU and GPU)
  /// @param ns       full size on the CPU
  /// @param streams  stream array for asynchronous transfer/computation
  /// @param ssizes   sizes of the slices to be put/get on GPU
  void resize(T init,
	      size_t ns, 
	      std::vector<cudaStream_t>* streams = NULL,
	      std::vector<size_t>* ssizes = NULL) 
  { 
    if (ns > 0) {
      cfassert(m_stride > 0);
      m_alloc.allocate(ns*m_stride, streams, ssizes);
      for (size_t i = 0; i < m_alloc.size(); ++i) {(*this)[i] = init;}
    }
  }
  
  /// resize storage (on CPU and GPU)
  /// @param ns       full size on the CPU
  void resize(size_t ns) 
  { 
    if (ns > 0) {
      cfassert(m_stride > 0);
      m_alloc.allocate(ns*m_stride, NULL, NULL);
    }
  }
  
  /// factor for which dividing the element size for compatibility with @see ParVector
  size_t sizeFactor() const {return sizeof(T);}
  
  /// get the stride
  size_t stride() const {return m_stride;}
  
  /// this is just kept momentarily for compatibility with ParVector 
  /// stride is the stride corresponding to the size of each element 
  void initialize (const T & Init, size_t initSize, size_t stride) {m_stride = stride;}
  
  /// synchronous copy of all data from CPU to GPU
  /// @post the device array will be allocated if not available yet
  void put() {CudaCopy<SELF, TO_GPU>(m_alloc.ptrDev(), m_alloc.ptr(), m_alloc.size());}
  
  /// synchronous copy of all data from GPU to CPU
  void get() {CudaCopy<SELF, FROM_GPU>(m_alloc.ptr(), m_alloc.ptrDev(), m_alloc.size());}
  
  /// get the raw pointer for accessing data on CPU
  T* ptr() {return m_alloc.ptr();}
  
  /// get the raw pointer for accessing data on GPU
  T* ptrDev() {return m_alloc.ptrDev();}
  
  /// subscripting operator for accessing data on CPU
  T operator[] (size_t i) const {cfassert(i < m_alloc.size()); return m_alloc.getv(i);}
  
  /// subscripting operator for accessing data on CPU
  T& operator[] (size_t i) {cfassert(i < m_alloc.size()); return m_alloc.getr(i);}
  
  /// operator for accessing data on CPU taking into account the stride
  /// this operator must be used judiciously
  T operator() (size_t i) const {cfassert(i < size()); return m_alloc.getv(i*m_stride);}
  
  /// operator for accessing data on CPU taking into account the stride
  /// this operator must be used judiciously
  T& operator() (size_t i) {cfassert(i < size()); return m_alloc.getr(i*m_stride);}
  
  /// overloading of the operator() returns a slice corresponding to the given stream ID
  CudaVectorSlice<T, ALLOC> operator() (int streamID)
  { 
    // here you need to allocate the device array
    return CudaVectorSlice<T,ALLOC>(m_alloc.size(streamID),
				    m_alloc.ptr(streamID), 
				    m_alloc.ptrDevS(streamID), 
				    m_alloc.stream(streamID));
  }
  
  /// size of the whole storage on the CPU
  /// @post WATCH OUT: this size != m_alloc.size() when m_stride != 1 (GLOBAL arrays) 
  size_t size() const {return m_alloc.size()/m_stride;}
  
  /// size of the slice corresponding to the given streamID
  size_t size(size_t streamID) const {return m_alloc.size(streamID);}
  
  /// size of the slice corresponding to the given streamID
  size_t start(size_t streamID) const {return m_alloc.start(streamID);}
  
private:
  
  /// allocator
  ALLOC<T> m_alloc;
  
  /// allocator
  size_t m_stride;
};

//////////////////////////////////////////////////////////////////////////////

/// This class implements a CUDA memory copy functionality via a computational 
/// constructor suitable for the two given template arguments. By default it 
/// does nothing.
///
/// @author Andrea Lani
template <typename ARRAY, TransferType TT> class CudaCopy {};

//////////////////////////////////////////////////////////////////////////////

/// This partial specialization of CudaCopy copies memory from host to device
/// synchronously or asynchronously depending on the given cudaStream_t pointer.
///
/// @author Andrea Lani
template <typename T> 
class CudaCopy<CudaVectorSlice<T, PinnedHostAlloc>, TO_GPU > {
public:
  CudaCopy(T* out, T* in, size_t ns, cudaStream_t* cs) 
  {
    (cs != NULL) ? CudaEnv::copyAsyncHost2Dev(out, in, ns, *cs) : 
      CudaEnv::copyHost2Dev(out, in, ns);  
  }
};

//////////////////////////////////////////////////////////////////////////////

/// This partial specialization of CudaCopy copies memory from device to host
/// synchronously or asynchronously depending on the given cudaStream_t pointer.
///
/// @author Andrea Lani
template <typename T> 
class CudaCopy<CudaVectorSlice<T, PinnedHostAlloc>, FROM_GPU > {
public:
  CudaCopy(T* out, T* in, size_t ns, cudaStream_t* cs) 
  {
    (cs != NULL) ? CudaEnv::copyAsyncDev2Host(out, in, ns, *cs) : 
      CudaEnv::copyDev2Host(out, in, ns); 
  }
};

//////////////////////////////////////////////////////////////////////////////

/// This partial specialization of CudaCopy copies memory from host to device
/// synchronously.
///
/// @author Andrea Lani
template <typename ARRAY> 
class CudaCopy<ARRAY, TO_GPU> {
public:
  CudaCopy(typename ARRAY::TYPE* out, typename ARRAY::TYPE* in, size_t ns, cudaStream_t* cs = NULL) 
  {
    cfassert(cs == NULL);
    CudaEnv::copyHost2Dev(out, in, ns);  
  }
};

//////////////////////////////////////////////////////////////////////////////

/// This partial specialization of CudaCopy copies memory from device to host
/// synchronously.
///
/// @author Andrea Lani
template <typename ARRAY> 
class CudaCopy<ARRAY, FROM_GPU> {
public:
  CudaCopy(typename ARRAY::TYPE* out, typename ARRAY::TYPE* in, size_t ns, cudaStream_t* cs = NULL) 
  {
    cfassert(cs == NULL);
    CudaEnv::copyDev2Host(out, in, ns);  
  }
};

//////////////////////////////////////////////////////////////////////////////
  
  } // namespace CudaEnv

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
