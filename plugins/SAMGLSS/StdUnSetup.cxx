// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SAMGLSS/StdUnSetup.hh"
#include "SAMGLSS/SAMGLSSModule.hh"
#include "Framework/MethodCommandProvider.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

Framework::MethodCommandProvider< StdUnSetup,SAMGLSSData,SAMGLSSModule >
  stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;
  getMethodData().getSolVector().destroy();
  getMethodData().getRhsVector().destroy();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SAMGLSS
} // namespace COOLFluiD

