// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "CrankNichSetup.hh"
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

MethodCommandProvider<CrankNichSetup, NewtonIteratorData, NewtonMethodModule> crankNichSetupProvider("CrankNichSetup");

//////////////////////////////////////////////////////////////////////////////

CrankNichSetup::CrankNichSetup(std::string name) : StdSetup(name),
  socket_pastRhs("pastRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

CrankNichSetup::~CrankNichSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
CrankNichSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_pastRhs);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CrankNichSetup::execute()
{
  CFAUTOTRACE;

  StdSetup::execute();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle<CFreal> pastRhs = socket_pastRhs.getDataHandle();
  pastRhs.resize(rhs.size());
  pastRhs = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
