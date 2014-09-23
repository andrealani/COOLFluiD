// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "BackwardEuler/BackwardEuler.hh"
#include "StdSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, BwdEulerData, BackwardEulerModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(std::string name) :
  BwdEulerCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbStates*nbEqs);
  rhs = 0.0;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();
  cf_assert(nbLSS > 0);
  updateCoeff.resize(nbStates*nbLSS);

  DataHandle<std::valarray<State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();
  bStatesNeighbors.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD
