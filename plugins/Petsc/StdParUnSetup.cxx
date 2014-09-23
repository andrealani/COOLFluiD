// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"

#include "Petsc/Petsc.hh"
#include "Petsc/StdParUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdParUnSetup, PetscLSSData, PetscModule> stdParUnSetupProvider("StdParUnSetup");

//////////////////////////////////////////////////////////////////////////////

void StdParUnSetup::execute()
{
  CFAUTOTRACE;

  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();
  KSP& ksp = getMethodData().getKSP();

  CF_CHKERRCONTINUE(KSPDestroy(&ksp));

  sol.destroy();
  rhs.destroy();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
