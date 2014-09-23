// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "MathTools/MathConsts.hh"
#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"
#include "SimpleGlobalMeshAdapter/DummyMeshInterpolator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DummyMeshInterpolator, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> DummyMeshInterpolatorProvider("DummyMeshInterpolator");

//////////////////////////////////////////////////////////////////////////////

DummyMeshInterpolator::DummyMeshInterpolator(const std::string& name)  :
  SimpleMeshAdapterCom(name),
  socket_states("states"),
  socket_otherStates("states")
{
}

//////////////////////////////////////////////////////////////////////////////

void DummyMeshInterpolator::configure ( Config::ConfigArgs& args )
{
  SimpleMeshAdapterCom::configure(args);

  socket_otherStates.setDataSocketNamespace(getMethodData().getOtherNamespace());
}

//////////////////////////////////////////////////////////////////////////////

void DummyMeshInterpolator::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbOtherStates = otherStates.size();

  RealVector result(PhysicalModelStack::getActive()->getNbEq());
  RealVector vector(PhysicalModelStack::getActive()->getDim());
  for(CFuint iState=0; iState< nbStates; iState++)
  {
    RealVector stateCoord = states[iState]->getCoordinates();
    CFreal minDistance = MathTools::MathConsts::CFrealMax();
    for(CFuint iOtherState=0; iOtherState< nbOtherStates; iOtherState++)
    {
      vector = stateCoord - (otherStates[iOtherState])->getCoordinates();
      CFreal distance = vector.norm2();

      if(distance < minDistance){
        minDistance = distance;
//         result = *(otherStates[iOtherState]);
        (states[iState])->copyData(*(otherStates[iOtherState]));
      }
    }
  }

//   for(CFuint iState=0; iState< nbStates; iState++)
//   {
//     *(states[iState]) = 2.;
//   }

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  DummyMeshInterpolator::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_states);
  result.push_back(&socket_otherStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
