// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Paralution/StdUnSetup.hh"
#include "Paralution/Paralution.hh"
#include "Framework/MethodCommandProvider.hh"

namespace COOLFluiD {
  namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider< StdUnSetup,ParalutionLSSData,ParalutionModule >
  stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  // destroy system vectors
  getMethodData().getSolVector().destroy();
  getMethodData().getRhsVector().destroy();

  // termination and release of memory
  getMethodData().unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Paralution
} // namespace COOLFluiD

