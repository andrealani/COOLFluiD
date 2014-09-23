// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ParaViewWriter/ParaViewWriter.hh"
#include "StdUnSetup.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  using namespace Framework;

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup, ParaWriterData, ParaViewWriterModule> stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD
