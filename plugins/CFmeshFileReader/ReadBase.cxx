// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <set>

#include "MathTools/MathChecks.hh"

#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"

#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/ReadBase.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

ReadBase::ReadBase(const std::string& name) :
  CFmeshReaderCom(name),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states"),
  m_data(new CFmeshReaderSource())
{
}

//////////////////////////////////////////////////////////////////////////////

ReadBase::~ReadBase()
{
}

//////////////////////////////////////////////////////////////////////////////

void ReadBase::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > > extraNVars =
    getMethodData().getExtraNVarSocketNamesAndTags();
  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > > extraSVars =
    getMethodData().getExtraSVarSocketNamesAndTags();
  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > > extraVars =
    getMethodData().getExtraVarSocketNamesAndTags();

  for(CFuint iNVar = 0; iNVar < extraNVars->size();iNVar++){
    std::string socketName = (*extraNVars)[iNVar].first;
    m_sockets.createSocketSource<CFreal>(socketName);
  }

  for(CFuint iSVar = 0; iSVar < extraSVars->size();iSVar++){
    std::string socketName = (*extraSVars)[iSVar].first;
    m_sockets.createSocketSource<CFreal>(socketName);
  }
  
  for(CFuint iVar = 0; iVar < extraVars->size();iVar++){
    std::string socketName = (*extraVars)[iVar].first;
    m_sockets.createSocketSource<CFreal>(socketName);
  }

  bool readPastNodes = getMethodData().readPastNodes();
  bool readPastStates = getMethodData().readPastStates();

  if(readPastNodes) m_sockets.createSocketSink<Node*>("pastNodes");
  if(readPastStates) m_sockets.createSocketSink<State*>("pastStates");

  
  bool readInterNodes = getMethodData().readInterNodes();
  bool readInterStates = getMethodData().readInterStates();

  if(readInterNodes) m_sockets.createSocketSink<Node*>("interNodes");
  if(readInterStates) m_sockets.createSocketSink<State*>("interStates");

}

//////////////////////////////////////////////////////////////////////////////

void ReadBase::correctStates()
{
  string stateSolutionFile = getMethodData().getStateSolutionFile();
  
  if (stateSolutionFile != "Null") {
    const CFuint rank = PE::GetPE().GetRank("Default");
    string infile = stateSolutionFile + ".P" + StringOps::to_str(rank);
    
    ifstream fin(infile.c_str());
    CFuint nbStates = 0;
    CFuint nbEqs = 0;
    fin >> nbStates >> nbEqs;
  
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
    cf_always_assert(nbStates == states.size());
    cf_always_assert(nbEqs == PhysicalModelStack::getActive()->getNbEq());

    CFuint globalID = 0;
    for (CFuint iState = 0; iState < states.size(); ++iState)
    {
      fin >> globalID;
      cf_always_assert(states[iState]->getGlobalID() == globalID);
      fin >> *states[iState];
     }  
  } 
  
  if(m_data->getNbEquations() != PhysicalModelStack::getActive()->getNbEq())
  {
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();

    for (CFuint iState = 0; iState < states.size(); ++iState)
    {
      states[iState]->resize(PhysicalModelStack::getActive()->getNbEq());
      *states[iState] = 0.0;
    }

    m_data->setNbEquations(PhysicalModelStack::getActive()->getNbEq());
  }

  const vector<CFreal>& sv = getMethodData().getStateScalingValues();
  if (sv.size() > 0) {
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    cf_always_assert(sv.size() == nbEqs);
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();

    for (CFuint iState = 0; iState < states.size(); ++iState) {
      for (CFuint i = 0; i < nbEqs; ++i) {
	(*states[iState])[i] /= sv[i];
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

void ReadBase::applyScalings()
{
  DataHandle<Node*, GLOBAL> nodes = socket_nodes.getDataHandle();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  std::vector<CFreal> scalingVector(dim);

  // this is used to rescale the given mesh and used it rescaled from now on
  if(getMethodData().isAnisotropicScaling())
  {
    scalingVector = getMethodData().getAnisotropicScalingVector();
  }
  else
  {
    for(CFuint iDim = 0; iDim < dim; ++iDim)
    {
      scalingVector[iDim] = getMethodData().getMeshScalingFactor();
    }
  }

  // if reference lenght is to close to zero, throw an exception
  for(CFuint iDim = 0; iDim < dim; ++iDim){
    if ( MathTools::MathChecks::isZero(scalingVector[iDim]) )
      throw Common::BadValueException (FromHere(),"Scaling factor is too close to zero");
  }

  // only apply it if the reference lenght is different enough from 1.0
  for(CFuint iDim = 0; iDim < dim; ++iDim){
    if ( MathTools::MathChecks::isNotEqual(scalingVector[iDim],1.0) )
    {
      for (CFuint i = 0; i < nodes.size(); ++i)
      {
        (*nodes[i])[iDim] /= scalingVector[iDim];
      }
    }
  }

  // Nodes are here adimensionalized with a reference length
  // this is used to rescale the given mesh but only during the computation,
  // to make adimensional computations, but it must be rescaled back to the original
  // size before I/O operations
  const CFreal refLength =
    PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // if reference lenght is to close to zero, throw an exception
  if ( MathTools::MathChecks::isZero(refLength) )
    throw Common::BadValueException (FromHere(),"Reference lenght is too close to zero");

  // only apply it if the reference lenght is different enough from 1.0
  if ( MathTools::MathChecks::isNotEqual(refLength,1.0) )
  {
    for (CFuint i = 0; i < nodes.size(); ++i)
    {
      *nodes[i] /= refLength;
    }
  }

  vector<CFuint> nodeSwitchIDs = getMethodData().getNodeSwitchIDs();
  RealVector tmpNodeBkp(dim);
  
  if (nodeSwitchIDs.size() == dim) {
      for (CFuint i = 0; i < nodes.size(); ++i)
      {
	  RealVector& node = *nodes[i];
	  tmpNodeBkp = node;
	  for (CFuint iDim = 0; iDim < dim; ++iDim)
	  {
	      node[iDim] = tmpNodeBkp[nodeSwitchIDs[iDim]];
	  }
      }
  }
  
  if(getMethodData().readPastNodes())
  {
    DataHandle<Node*> pastNodes = m_sockets.getSocketSink<Node*>("pastNodes")->getDataHandle();
    //Apply Scalings to the pastNodes also
    for(CFuint iDim = 0; iDim < dim; ++iDim){
      if ( MathTools::MathChecks::isNotEqual(scalingVector[iDim],1.0) )
      {
        for (CFuint i = 0; i < nodes.size(); ++i)
        {
          (*pastNodes[i])[iDim] /= scalingVector[iDim];
        }
      }
    }

    // only apply it if the reference lenght is different enough from 1.0
    if ( MathTools::MathChecks::isNotEqual(refLength,1.0) )
    {
      for (CFuint i = 0; i < nodes.size(); ++i)
      {
        *pastNodes[i] /= refLength;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

void ReadBase::applyTranslation()
{
  DataHandle<Node*,GLOBAL> nodes = socket_nodes.getDataHandle();

  if (getMethodData().isTranslated())
  {
    std::vector<CFreal> translationVector = getMethodData().getTranslationVector();

    cf_assert(translationVector.size() == (*nodes[0]).size());
    for (CFuint i = 0; i < nodes.size(); ++i)
    {
      for (CFuint j = 0; j < translationVector.size(); ++j)
      {
        (*nodes[i])[j] += translationVector[j];
      }
    }

    if(getMethodData().readPastNodes())
    {
      DataHandle<Node*> pastNodes = m_sockets.getSocketSink<Node*>("pastNodes")->getDataHandle();
      cf_assert(translationVector.size() == (*pastNodes[0]).size());
      for (CFuint i = 0; i < pastNodes.size(); ++i)
      {
        for (CFuint j = 0; j < translationVector.size(); ++j)
        {
          (*pastNodes[i])[j] += translationVector[j];
        }
      }

    }

  }



}

/////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ReadBase::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

/////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ReadBase::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD
