// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Petsc_PetscHeaders_hh
#define COOLFluiD_Petsc_PetscHeaders_hh

//////////////////////////////////////////////////////////////////////////////

// must come first because petsc gets confused with some definitions of boost
#include <boost/filesystem/path.hpp>
#include "Common/PE.hh"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscerror.h"

// undefine the restrict as defined by petsc
#undef restrict

#ifndef NDEBUG
#ifndef CF_HAVE_CUDA // this is for petsc-dev
#define CF_CHKERRCONTINUE(n) if (n) {PetscError(Common::PE::GetPE().GetCommunicator(),__LINE__,__FUNCT__,__FILE__,__SDIR__,n,(PetscErrorType)0," ");}
#else
#define CF_CHKERRCONTINUE(n) if (n) {PetscError(Common::PE::GetPE().GetCommunicator(),__LINE__,__FUNCT__,__FILE__,n,(PetscErrorType)0," ");}
#endif
#else
#define CF_CHKERRCONTINUE(n) { n; }
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Petsc_PetscMatrix_hh
