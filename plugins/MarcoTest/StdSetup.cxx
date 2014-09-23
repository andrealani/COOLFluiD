// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MarcoTest/MarcoTest.hh"
#include "MarcoTest/MarcoTestMethodData.hh"
#include "Framework/BaseSetupFVMCC.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////
      
MethodCommandProvider<BaseSetupFVMCC<MarcoTestMethodCom>,
                      MarcoTestMethodData,
                      MarcoTestModule>
stdSetupMarcoTestProvider("StdSetup");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest

} // namespace COOLFluiD

