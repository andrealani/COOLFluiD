// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/StdUnSetup.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  using namespace Framework;

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup, TecWriterData, TecplotWriterModule> stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(std::string name) : TecWriterCom(name) {}

StdUnSetup::~StdUnSetup() {}

void StdUnSetup::execute() { CFAUTOTRACE; }

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
