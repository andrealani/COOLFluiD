// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/StdSetup.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,FluxReconstructionSolverData,FluxReconstructionModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states"),
  socket_gradients("gradients"),
  socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  StdSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_gradients);
  result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "StdSetup::execute() => START\n");

  // DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  // DataHandle<ProxyDofIterator< RealVector >* > nstatesProxy =
  //   socket_nstatesProxy.getDataHandle();

  // const CFuint nbStates = states.size();
  // nstatesProxy.resize(1);

  // // set node to state mapping
  // m_nodeIDToStateID.resize(nbStates);
  // for (CFuint stateID=0;stateID<nbStates;++stateID) {
  //   const CFuint nodeID = states[stateID]->getCoordinates().getLocalID();
  //   cf_assert(nodeID < nbStates);
  //   m_nodeIDToStateID[nodeID] = stateID;
  // }
  // nstatesProxy[0] =
  //   new DofDataHandleIterator< RealVector,State, GLOBAL >(states,&m_nodeIDToStateID);
  
    
  //   ///
  // // number of equations
  // const CFuint nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // // get number of gradients (assumed equal to the number of physical variables)
  // const CFuint nbrGrads = nbrEqs;
  
  // // dimensionality
  // const CFuint dim = PhysicalModelStack::getActive()->getDim();  

  // // get datahandle
  // DataHandle< std::vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // // resize gradients
  // gradients.resize(nbStates);
  // for (CFuint iState = 0; iState < nbStates; ++iState)
  // {
  //   gradients[iState].resize(nbrGrads);
  //   for (CFuint iGrad = 0; iGrad < nbrGrads; ++iGrad)
  //   {
  //     gradients[iState][iGrad].resize(dim);
  //   }
  // }
  /// 
  
  CFLog(VERBOSE, "StdSetup::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

