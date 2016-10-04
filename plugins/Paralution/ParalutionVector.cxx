// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Framework/MeshData.hh"
#include "Paralution/ParalutionVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

ParalutionVector::ParalutionVector() :
  Framework::LSSVector(),
  m_toBeDestroyed(false)
{
}

//////////////////////////////////////////////////////////////////////////////

ParalutionVector::~ParalutionVector()
{
  destroy();
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void ParalutionVector::create(MPI_Comm comm,
			 const CFint m,
			 const CFint M,
			 const char* name)
{  
  // // creation of the solution vector
//   CF_CHKERRCONTINUE( VecCreate(comm, &m_vec) );
  
// #ifdef CF_HAVE_CUDA
//   if (m_useGPU) {
//     CF_CHKERRCONTINUE( VecSetType(m_vec, VECCUSP) ); 
//   }
//   else
// #endif
//     {CF_CHKERRCONTINUE( VecSetType(m_vec, VECMPI) );}
  
//   CF_CHKERRCONTINUE( VecSetSizes(m_vec, m, M) );
  
//   CF_CHKERRCONTINUE( VecSetFromOptions(m_vec) );

//   CF_CHKERRCONTINUE( ParalutionObjectSetName((ParalutionObject) m_vec, name) );

  m_toBeDestroyed = true;
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void ParalutionVector::initialize(MPI_Comm comm,
				  const CFreal value)
{
  cf_assert(m_toBeDestroyed == true);
}
#endif
      
//////////////////////////////////////////////////////////////////////////////

void ParalutionVector::destroy()
{
  if (m_toBeDestroyed)
  {
    // CF_CHKERRCONTINUE( VecDestroy(&m_vec) );
  }

  m_toBeDestroyed = false;
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionVector::printToScreen() const
{

  // CF_CHKERRCONTINUE( VecView(m_vec, PETSC_VIEWER_STDOUT_WORLD) );

}

//////////////////////////////////////////////////////////////////////////////

void ParalutionVector::printToFile(const char* fileName) const
{
  // const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  
  // ParalutionViewer viewer;
  // CF_CHKERRCONTINUE( ParalutionViewerASCIIOpen(Common::PE::GetPE().GetCommunicator(nsp),
  // 					  fileName,
  // 					  &viewer) );
  // CF_CHKERRCONTINUE( VecView(m_vec, viewer) );
  
  // CF_CHKERRCONTINUE( ParalutionViewerDestroy(&viewer) );
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionVector::getValues(const CFuint m, const CFint* im, CFreal* values)
{
  // CFreal* vecValues;
  // CF_CHKERRCONTINUE( VecGetArray(m_vec, &vecValues) );

  // for (CFuint i=0; i<m; ++i) {
  //   values[i] = vecValues[im[i]];
  // }
  // CF_CHKERRCONTINUE( VecRestoreArray(m_vec, &vecValues) );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
