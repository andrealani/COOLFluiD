// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Common/PE.hh"
#include "Framework/MeshData.hh"
#include "Petsc/PetscVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

PetscVector::PetscVector() :
  Framework::LSSVector(),
  m_vec(),
  m_toBeDestroyed(false)
{
}

//////////////////////////////////////////////////////////////////////////////

PetscVector::~PetscVector()
{
  destroy();
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PetscVector::create(MPI_Comm comm,
			 const CFint m,
			 const CFint M,
			 const char* name)
{  
  // creation of the solution vector
  CF_CHKERRCONTINUE( VecCreate(comm, &m_vec) );
  
#ifdef CF_HAVE_CUDA
  if (m_useGPU) {
#ifdef CF_HAVE_VIENNACL
    CF_CHKERRCONTINUE( VecSetType(m_vec, VECVIENNACL) );
#else
#if PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 
    CF_CHKERRCONTINUE( VecSetType(m_vec, VECCUDA) ); 
#else
    CF_CHKERRCONTINUE( VecSetType(m_vec, VECCUSP) ); 
#endif
#endif 
  }
  else
#endif
    {CF_CHKERRCONTINUE( VecSetType(m_vec, VECMPI) );}
  
  CF_CHKERRCONTINUE( VecSetSizes(m_vec, m, M) );
  
  CF_CHKERRCONTINUE( VecSetFromOptions(m_vec) );
  
  CF_CHKERRCONTINUE( PetscObjectSetName((PetscObject) m_vec, name) );
  
  m_toBeDestroyed = true;
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PetscVector::initialize(MPI_Comm comm,
			     const CFreal value)
{
  cf_assert(m_toBeDestroyed == true);

  // check if random initialization is chosen
  PetscBool flg;
  
  // AL: this needs to be checked for PETSC 3.6.3
#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  flg = PETSC_FALSE;
#else
  CF_CHKERRCONTINUE( PetscOptionsHasName(PETSC_NULL,"-random_exact_sol",&flg) );
#endif
  
  if (flg) {
    PetscRandom rctx;
    CF_CHKERRCONTINUE( PetscRandomCreate(comm, &rctx) );

    CF_CHKERRCONTINUE( VecSetRandom(m_vec, rctx) );

    CF_CHKERRCONTINUE( PetscRandomDestroy(&rctx) );
  }
  else {
    CF_CHKERRCONTINUE( VecSet(m_vec, value) );
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////

void PetscVector::duplicate(PetscVector& other) const
{
  cf_assert(m_toBeDestroyed == true);

  CF_CHKERRCONTINUE( VecDuplicate(m_vec, &other.m_vec) );

  other.m_toBeDestroyed = true;
}

//////////////////////////////////////////////////////////////////////////////

void PetscVector::destroy()
{
  if (m_toBeDestroyed)
  {
    CF_CHKERRCONTINUE( VecDestroy(&m_vec) );
  }

  m_toBeDestroyed = false;
}

//////////////////////////////////////////////////////////////////////////////

void PetscVector::printToScreen() const
{

  CF_CHKERRCONTINUE( VecView(m_vec, PETSC_VIEWER_STDOUT_WORLD) );

}

//////////////////////////////////////////////////////////////////////////////

void PetscVector::printToFile(const char* fileName) const
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  
  PetscViewer viewer;
  CF_CHKERRCONTINUE( PetscViewerASCIIOpen(Common::PE::GetPE().GetCommunicator(nsp),
					  fileName,
					  &viewer) );
  CF_CHKERRCONTINUE( VecView(m_vec, viewer) );
  
  CF_CHKERRCONTINUE( PetscViewerDestroy(&viewer) );
}

//////////////////////////////////////////////////////////////////////////////

void PetscVector::getValues(const CFuint m, const CFint* im, CFreal* values)
{
  CFreal* vecValues;
  CF_CHKERRCONTINUE( VecGetArray(m_vec, &vecValues) );

  for (CFuint i=0; i<m; ++i) {
    values[i] = vecValues[im[i]];
  }
  CF_CHKERRCONTINUE( VecRestoreArray(m_vec, &vecValues) );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
