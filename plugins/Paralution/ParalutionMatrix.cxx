// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Framework/BlockAccumulator.hh"
#include "Paralution/ParalutionMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

ParalutionMatrix::ParalutionMatrix() :
  Framework::LSSMatrix()
{
}
      
//////////////////////////////////////////////////////////////////////////////

ParalutionMatrix::~ParalutionMatrix()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::createSeqAIJ(const CFint m,
                               const CFint n,
                               const CFint nz,
                               const CFint* nnz,
                               const char* name)
{
  // // creation of the matrix
  // CF_CHKERRCONTINUE(MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, nz, nnz, &m_mat));

  // if (name != CFNULL) {
  //   CF_CHKERRCONTINUE(ParalutionObjectSetName((ParalutionObject) m_mat, name));
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::createSeqBAIJ(const CFuint blockSize,
                                const CFint m,
                                const CFint n,
                                const CFint nz,
                                const CFint* nnz,
                                const char* name)
{
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void ParalutionMatrix::createParAIJ(MPI_Comm comm,
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
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void ParalutionMatrix::createParBAIJ(MPI_Comm comm,
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
}
#endif

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::setValues(const Framework::BlockAccumulator& acc)
{
  // CFLog(DEBUG_MIN, "ParalutionMatrix::setValues()\n");
  // CF_CHKERRCONTINUE( MatSetValuesBlocked(m_mat,acc.getM(),&acc.getIM()[0],acc.getN(),&acc.getIN()[0],
  //   const_cast<Framework::BlockAccumulator&>(acc).getPtr(), INSERT_VALUES) );
}
      
//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::addValues(const Framework::BlockAccumulator& acc)
{
  // CFLog(DEBUG_MIN, "ParalutionMatrix::addValues()\n");
  // CF_CHKERRCONTINUE( MatSetValuesBlocked(m_mat,acc.getM(),&acc.getIM()[0],acc.getN(),&acc.getIN()[0],
  // 					 const_cast<Framework::BlockAccumulator&>(acc).getPtr(), ADD_VALUES) );
}
      
//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::printToScreen() const
{
  // CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY));
  // CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY));
  // CF_CHKERRCONTINUE(MatView(m_mat, PETSC_VIEWER_STDOUT_WORLD));
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::printToFile(const char* fileName) const
{
  // ParalutionViewer viewer;
  // CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY));
  // CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY));
  
  // const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  // CF_CHKERRCONTINUE(ParalutionViewerASCIIOpen(Common::PE::GetPE().GetCommunicator(nsp),fileName,&viewer));
  // CF_CHKERRCONTINUE(ParalutionViewerSetType(viewer,PETSCVIEWERASCII));
  // CF_CHKERRCONTINUE(ParalutionViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATRIXMARKET));
  // // PETSC_VIEWER_ASCII_COMMON, PETSC_VIEWER_ASCII_INDEX, PETSC_VIEWER_ASCII_MATLAB
  // CF_CHKERRCONTINUE(MatView(m_mat, viewer));
  // CF_CHKERRCONTINUE(ParalutionViewerDestroy(&viewer));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
