// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTypes_hh
#define COOLFluiD_MathTypes_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CFVec.hh"
#include "Common/CUDA/CFVecSlice.hh"
#include "Common/CUDA/CFMat.hh"
#include "Common/CUDA/CFMatSlice.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////
  
/**
 * This class defines mathematical types to be used in algorithmic kernels.
 * By default CPU devices are assumed.
 *
 * @author Andrea Lani
 *
 */
template <typename T, DeviceType DT = CPU, int SIZE = 0> 
class MathTypes
{
public:
  enum {N=SIZE};
  typedef MathTools::CFVec<T,SIZE> VEC;
  typedef MathTools::CFMat<T,SIZE,SIZE> MAT;
  typedef MathTools::CFVecSlice<T,SIZE> SLICEVEC;
  typedef MathTools::CFMatSlice<T,SIZE,SIZE> SLICEMAT;
}; // class MathTypes
  
#ifdef CF_HAVE_CUDA
/**
 * This class defines mathematical types to be used in algorithmic kernels.
 * This is a partial specialization for GPU devices.
 *
 * @author Andrea Lani
 *
 */
template <typename T, int SIZE> 
class MathTypes<T, GPU, SIZE>
{
public:
  enum {N=SIZE};
  typedef CudaEnv::CFVec<T,SIZE> VEC;
  typedef CudaEnv::CFMat<T,SIZE,SIZE> MAT;
  typedef CudaEnv::CFVecSlice<T,SIZE> SLICEVEC;
  typedef CudaEnv::CFMatSlice<T,SIZE,SIZE> SLICEMAT;
}; // class MathTypes
#endif
  
///////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTypes_hh
