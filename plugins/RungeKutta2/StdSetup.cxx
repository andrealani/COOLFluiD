// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKutta2/RungeKutta2.hh"
#include "StdSetup.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, RK2Data, RungeKutta2Module> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  RK2Com(name),
  socket_rhs("rhs"),
  socket_u0("u0"),
  socket_updateCoeff("updateCoeff"),
  socket_upCoeff0("upCoeff0"),
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
  result.push_back(&socket_u0);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_upCoeff0);

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

  DataHandle<RealVector> u0 = socket_u0.getDataHandle();
  u0.resize(nbStates);
  for (CFuint i = 0; i < nbStates; i++) {
    u0[i].resize(nbEqs);
  }

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(nbStates);

  DataHandle<CFreal> upCoeff0 = socket_upCoeff0.getDataHandle();
  upCoeff0.resize(nbStates);
 }

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD
