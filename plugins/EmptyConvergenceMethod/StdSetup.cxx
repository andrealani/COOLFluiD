// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"

#include "EmptyConvergenceMethod/EmptyConvergenceMethod.hh"
#include "EmptyConvergenceMethod/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace EmptyConvergenceMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, EmptyIteratorData, EmptyConvergenceMethodLib> 
stdSetupEmptyIteratorProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  EmptyIteratorCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::configure ( Config::ConfigArgs& args )
{
  EmptyIteratorCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;
  
  DataHandle < State*, GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbStates*nbEqs);
  rhs = 0.0;
  
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace EmptyConvergenceMethod

} // namespace COOLFluiD
