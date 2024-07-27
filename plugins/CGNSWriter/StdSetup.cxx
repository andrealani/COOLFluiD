// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"

#include "CGNSWriter/CGNSWriter.hh"
#include "CGNSWriter/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, CGWriterData, CGNSWriterModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(std::string name) : CGWriterCom(name) {}

StdSetup::~StdSetup() {}

void StdSetup::execute() { CFAUTOTRACE; }

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD
