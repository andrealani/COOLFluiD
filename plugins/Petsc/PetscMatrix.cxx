// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Common/PE.hh"
#include "Framework/BlockAccumulator.hh"
#include "Petsc/PetscMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

PetscMatrix::PetscMatrix() :
  Framework::LSSMatrix(),
  m_mat(),
  _isMatShell(false),
  _isAIJ(false)
{
}
      
//////////////////////////////////////////////////////////////////////////////

PetscMatrix::~PetscMatrix()
{
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::createSeqAIJ(const CFint m,
                               const CFint n,
                               const CFint nz,
                               const CFint* nnz,
                               const char* name)
{
  // creation of the matrix
  CF_CHKERRCONTINUE(MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, nz, nnz, &m_mat));

  if (name != CFNULL) {
    CF_CHKERRCONTINUE(PetscObjectSetName((PetscObject) m_mat, name));
  }
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::createSeqBAIJ(const CFuint blockSize,
                                const CFint m,
                                const CFint n,
                                const CFint nz,
                                const CFint* nnz,
                                const char* name)
{
  CF_CHKERRCONTINUE(MatCreate(PETSC_COMM_SELF, &m_mat));
  CF_CHKERRCONTINUE(MatSetSizes(m_mat, m, n, m, n));
  
  if (_isAIJ) {
#ifdef CF_HAVE_CUDA
    if (m_useGPU) {
#ifdef CF_HAVE_VIENNACL
     CF_CHKERRCONTINUE(MatSetType(m_mat, MATSEQAIJVIENNACL));
#else
#if PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12
     CF_CHKERRCONTINUE(MatSetType(m_mat, MATSEQAIJCUSPARSE));
#else
     CF_CHKERRCONTINUE(MatSetType(m_mat, MATSEQAIJCUSP));
#endif
#endif
    }
    else 
#endif   
      {CF_CHKERRCONTINUE(MatSetType(m_mat, MATSEQAIJ));}

    CFLog(VERBOSE, "PetscMatrix::createSeqAIJ() on GPU["<< m_useGPU << "]\n" );   
 
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12
    CF_CHKERRCONTINUE(MatSetBlockSize(m_mat, blockSize));
    CF_CHKERRCONTINUE(MatSeqAIJSetPreallocation(m_mat, nz, nnz));
    CFLog(DEBUG_MIN, "PetscMatrix::createSeqBAIJ() => MatSeqAIJSetPreallocation() \n" );
#else
    CF_CHKERRCONTINUE(MatSeqAIJSetPreallocation(m_mat, nz, nnz)); 
    CF_CHKERRCONTINUE(MatSetBlockSize(m_mat, blockSize));
#endif
  }
  else {
    CF_CHKERRCONTINUE(MatSetType(m_mat, MATSEQBAIJ));
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12
    CF_CHKERRCONTINUE(MatSetBlockSize(m_mat, blockSize));
#endif
    CF_CHKERRCONTINUE(MatSeqBAIJSetPreallocation(m_mat, blockSize, nz, nnz));
  }
    
  if (name != CFNULL) {
    CF_CHKERRCONTINUE(PetscObjectSetName((PetscObject) m_mat, name));
  }
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::createSeqJFMat(const CFint m,
				 const CFint n,
				 const CFint M,
				 const CFint N,
				 void* ctx,
				 const char* name)
{
  // creation of the matrix
  CF_CHKERRCONTINUE(MatCreateShell(PETSC_COMM_SELF, m, n, M, N, ctx, &m_mat));
  
  _isMatShell = false;
  
  if (name != CFNULL) {
    CF_CHKERRCONTINUE(PetscObjectSetName((PetscObject) m_mat, name));
  }
}

// //////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PetscMatrix::createParAIJ(MPI_Comm comm,
                               const CFint m,
                               const CFint n,
                               const CFint M,
                               const CFint N,
                               const CFint dnz,
                               const CFint* dnnz,
                               const CFint onz,
                               const CFint* onnz,
                               const char* name)
{
 #if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
   CF_CHKERRCONTINUE(MatCreateAIJ(comm, m, n, M, N, dnz, dnnz,  onz, onnz, &m_mat));
#else
  CF_CHKERRCONTINUE(MatCreateMPIAIJ(comm, m, n, M, N, dnz, dnnz,  onz, onnz, &m_mat));
#endif
  if (name != CFNULL) {
    CF_CHKERRCONTINUE(PetscObjectSetName((PetscObject) m_mat, name));
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PetscMatrix::createParBAIJ(MPI_Comm comm,
				const CFuint blockSize,
				const CFint m,
				const CFint n,
				const CFint M,
				const CFint N,
				const CFint dnz,
				const CFint* dnnz,
				const CFint onz,
				const CFint* onnz,
				const char* name)
{
  CF_CHKERRCONTINUE(MatCreate(comm, &m_mat));
  CF_CHKERRCONTINUE(MatSetSizes(m_mat, m, n, M, N));

  if (_isAIJ) {
#ifdef CF_HAVE_CUDA
    if (m_useGPU) {
#if PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
      CF_CHKERRCONTINUE(MatSetType(m_mat, MATMPIAIJCUSPARSE));
#else
      CF_CHKERRCONTINUE(MatSetType(m_mat, MATMPIAIJCUSP));
#endif
    }
    else
#endif
      {CF_CHKERRCONTINUE(MatSetType(m_mat, MATMPIAIJ));}

    CFLog(VERBOSE, "PetscMatrix::createParAIJ() on GPU["<< m_useGPU << "]\n" );
    
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
    CF_CHKERRCONTINUE(MatSetBlockSize(m_mat, blockSize));
    CF_CHKERRCONTINUE(MatMPIAIJSetPreallocation(m_mat, dnz, dnnz, onz, onnz));
#else
    CF_CHKERRCONTINUE(MatMPIAIJSetPreallocation(m_mat, dnz, dnnz, onz, onnz));
    CF_CHKERRCONTINUE(MatSetBlockSize(m_mat, blockSize));
#endif
  }
  else {
    CF_CHKERRCONTINUE(MatSetType(m_mat, MATMPIBAIJ));
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
    CF_CHKERRCONTINUE(MatSetBlockSize(m_mat, blockSize));
#endif
    CF_CHKERRCONTINUE(MatMPIBAIJSetPreallocation(m_mat, blockSize, dnz, dnnz, onz, onnz));
  }
  
  if (name != CFNULL) {
    CF_CHKERRCONTINUE(PetscObjectSetName((PetscObject) m_mat, name));
  }
  
  // option to improve efficiency and reliability
  // CF_CHKERRCONTINUE(MatSetOption(m_mat, MAT_USE_HASH_TABLE, PETSC_TRUE)); 
  // this gives a wrong residual in parallel
  CF_CHKERRCONTINUE(MatSetOption(m_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE)); 
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PetscMatrix::createParJFMat(MPI_Comm comm,
				 const CFint m,
				 const CFint n,
				 const CFint M,
				 const CFint N,
				 void* ctx,
				 const char* name)
{
  // creation of the matrix
  CF_CHKERRCONTINUE(MatCreateShell(comm, m, n, M, N, ctx, &m_mat));
  _isMatShell = true;
  
  if (name != CFNULL) {
    CF_CHKERRCONTINUE(PetscObjectSetName((PetscObject) m_mat, name));
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::setValues(const Framework::BlockAccumulator& acc)
{
  CFLog(DEBUG_MIN, "PetscMatrix::setValues()\n");
  CF_CHKERRCONTINUE( MatSetValuesBlocked(m_mat,acc.getM(),&acc.getIM()[0],acc.getN(),&acc.getIN()[0],
    const_cast<Framework::BlockAccumulator&>(acc).getPtr(), INSERT_VALUES) );
}
      
//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::addValues(const Framework::BlockAccumulator& acc)
{
  CFLog(DEBUG_MIN, "PetscMatrix::addValues()\n");
  CF_CHKERRCONTINUE( MatSetValuesBlocked(m_mat,acc.getM(),&acc.getIM()[0],acc.getN(),&acc.getIN()[0],
					 const_cast<Framework::BlockAccumulator&>(acc).getPtr(), ADD_VALUES) );
}
      
//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::printToScreen() const
{
  CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY));
  CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY));
  CF_CHKERRCONTINUE(MatView(m_mat, PETSC_VIEWER_STDOUT_WORLD));
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::printToFile(const char* fileName) const
{
/*
  // PETSC VIEWER
  #define PETSC_VIEWER_SOCKET       "socket"
  #define PETSC_VIEWER_ASCII        "ascii"
  #define PETSC_VIEWER_BINARY       "binary"
  #define PETSC_VIEWER_STRING       "string"
  #define PETSC_VIEWER_DRAW         "draw"
  #define PETSC_VIEWER_VU           "vu"
  #define PETSC_VIEWER_MATHEMATICA  "mathematica"
  #define PETSC_VIEWER_SILO         "silo"
  #define PETSC_VIEWER_NETCDF       "netcdf"
  #define PETSC_VIEWER_HDF5         "hdf5"
  #define PETSC_VIEWER_MATLAB       "matlab"

  // PETSC VIEWER FORMAT
  PETSC_VIEWER_DEFAULT,
  PETSC_VIEWER_ASCII_MATLAB, 
  PETSC_VIEWER_ASCII_MATHEMATICA,
  PETSC_VIEWER_ASCII_IMPL,
  PETSC_VIEWER_ASCII_INFO,
  PETSC_VIEWER_ASCII_INFO_DETAIL,
  PETSC_VIEWER_ASCII_COMMON,
  PETSC_VIEWER_ASCII_SYMMODU,
  PETSC_VIEWER_ASCII_INDEX,
  PETSC_VIEWER_ASCII_DENSE,
  PETSC_VIEWER_ASCII_MATRIXMARKET,
  PETSC_VIEWER_ASCII_VTK,
  PETSC_VIEWER_ASCII_VTK_CELL,
  PETSC_VIEWER_ASCII_VTK_COORDS,
  PETSC_VIEWER_ASCII_PCICE,
  PETSC_VIEWER_ASCII_PYLITH,
  PETSC_VIEWER_ASCII_PYLITH_LOCAL,
  PETSC_VIEWER_ASCII_PYTHON,
  PETSC_VIEWER_ASCII_FACTOR_INFO,
  PETSC_VIEWER_DRAW_BASIC,
  PETSC_VIEWER_DRAW_LG,
  PETSC_VIEWER_DRAW_CONTOUR, 
  PETSC_VIEWER_DRAW_PORTS,
  PETSC_VIEWER_NATIVE,
  PETSC_VIEWER_NOFORMAT
*/
  PetscViewer viewer;
  CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY));
  CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY));
  
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  CF_CHKERRCONTINUE(PetscViewerASCIIOpen(Common::PE::GetPE().GetCommunicator(nsp),fileName,&viewer));
  CF_CHKERRCONTINUE(PetscViewerSetType(viewer,PETSCVIEWERASCII));
  CF_CHKERRCONTINUE(PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATRIXMARKET));
  // PETSC_VIEWER_ASCII_COMMON, PETSC_VIEWER_ASCII_INDEX, PETSC_VIEWER_ASCII_MATLAB
  CF_CHKERRCONTINUE(MatView(m_mat, viewer));
  CF_CHKERRCONTINUE(PetscViewerDestroy(&viewer));
  
/*
  CF_CHKERRCONTINUE(PetscViewerHDF5Open(PETSC_COMM_SELF,fileName,FILE_MODE_WRITE,&viewer));
  CF_CHKERRCONTINUE(MatView(m_mat, viewer));
  CF_CHKERRCONTINUE(PetscViewerDestroy(viewer));
*/
/*
  CF_CHKERRCONTINUE(PetscViewerMatlabOpen(PETSC_COMM_SELF,fileName,FILE_MODE_WRITE,&viewer));
  CF_CHKERRCONTINUE(MatView(m_mat, viewer));
  CF_CHKERRCONTINUE(PetscViewerDestroy(viewer));
*/
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::setJFFunction(void (*JFFunction)(void))
{
  CF_CHKERRCONTINUE(MatShellSetOperation(m_mat, MATOP_MULT, JFFunction));
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::getJFContext(void **ctx)
{
  CF_CHKERRCONTINUE(MatShellGetContext(m_mat, ctx));
}

//////////////////////////////////////////////////////////////////////////////

void PetscMatrix::setAIJ(bool isAIJ)
{
  _isAIJ = isAIJ;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
