// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CGNSWriter/CGNSWriter.hh"
#include "CGNSWriter/StdUnSetup.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  using namespace Framework;

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup, CGWriterData, CGNSWriterModule> stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(std::string name) : CGWriterCom(name) {}

StdUnSetup::~StdUnSetup() {}

void StdUnSetup::execute() { CFAUTOTRACE; }

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD
