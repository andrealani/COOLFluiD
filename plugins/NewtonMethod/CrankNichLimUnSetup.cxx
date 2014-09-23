// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "CrankNichLimUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CrankNichLimUnSetup, NewtonIteratorData, NewtonMethodModule> CrankNichLimUnSetupProvider("CrankNichLimUnSetup");

//////////////////////////////////////////////////////////////////////////////

CrankNichLimUnSetup::CrankNichLimUnSetup(std::string name) : CrankNichUnSetup(name),
socket_pastPastStates("pastPastStates")
{
}

//////////////////////////////////////////////////////////////////////////////

CrankNichLimUnSetup::~CrankNichLimUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CrankNichLimUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = CrankNichUnSetup::needsSockets();

  result.push_back(&socket_pastPastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CrankNichLimUnSetup::execute()
{
  CFAUTOTRACE;

  CrankNichUnSetup::execute();

  DataHandle<State*> pastPastStates = socket_pastPastStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastPastStates.size(); ++i) {
    delete pastPastStates[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
