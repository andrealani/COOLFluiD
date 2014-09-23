// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "BDF2Setup.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BDF2Setup, NewtonIteratorData, NewtonMethodModule> BDF2SetupProvider("BDF2Setup");

//////////////////////////////////////////////////////////////////////////////

BDF2Setup::BDF2Setup(std::string name) : StdSetup(name),
  socket_pastTimeRhs("pastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

BDF2Setup::~BDF2Setup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
BDF2Setup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BDF2Setup::execute()
{
  CFAUTOTRACE;

  StdSetup::execute();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
  pastTimeRhs.resize(rhs.size());
  pastTimeRhs = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
