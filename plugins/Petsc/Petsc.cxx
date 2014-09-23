// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Common/CFLog.hh"
#include "Environment/CFEnv.hh"
#include "Environment/CFEnvVars.hh"

#include "Petsc/Petsc.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

void PetscModule::initiate()
{
  if (!isInitialized())
  {
    std::pair<int,char**> args = Environment::CFEnv::getInstance().getVars()->InitArgs;
    PetscInitialize(&args.first, &args.second, (char*)0, 0);
    CFLogNotice("Initiated PETSc\n");
    m_init = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void PetscModule::terminate()
{
  if(isInitialized())
  {
    PetscFinalize();
    CFLogNotice("Terminated PETSc\n");
    m_init = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
