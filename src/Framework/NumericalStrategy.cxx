// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Config/BadMatchException.hh"

#include "Framework/NumericalStrategy.hh"
#include "Framework/BaseDataSocketSource.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

NumericalStrategy::NumericalStrategy(const std::string& name)
  : Common::OwnedObject(),
    ConfigObject(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NumericalStrategy::~NumericalStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void NumericalStrategy::configureNestedSockets( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure the source sockets
  vector<Common::SafePtr<BaseDataSocketSource> >::iterator srcItr;
  vector<Common::SafePtr<BaseDataSocketSource> > providedSockets = providesSockets();
  for (srcItr = providedSockets.begin(); srcItr != providedSockets.end(); ++srcItr) {
      cf_assert(srcItr->isNotNull());
      ConfigObject& obj = *(*srcItr);
      configureNested ( obj, args );
  }

  // configure the sink sockets
  vector<Common::SafePtr<BaseDataSocketSink> >::iterator snkItr;
  vector<Common::SafePtr<BaseDataSocketSink> > neededSockets = needsSockets();
  for (snkItr = neededSockets.begin(); snkItr != neededSockets.end(); ++snkItr) {
      cf_assert(snkItr->isNotNull());
      configureNested(*(*snkItr), args);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NumericalStrategy::allocateStrategySockets()
{
  CFAUTOTRACE;
    
  typedef vector<Common::SafePtr<BaseDataSocketSource> > VecSock;
  VecSock providedSockets = providesSockets();

  if ( providedSockets.empty() ) return;

  // loop on the MeshData's to allocate the sockets in the
  // the one that defines the namespace of the socket
  typedef std::vector<Common::SafePtr<MeshData> > MDLst;
  MDLst mds = MeshDataStack::getInstance().getAllEntries();
  MDLst::iterator meshData = mds.begin();
  for(; meshData != mds.end(); ++meshData)
  {
    CF_DEBUG_OBJ ((*meshData)->getPrimaryNamespace());
    CF_DEBUG_OBJ (providedSockets.size());

    VecSock::iterator srcItr = providedSockets.begin();
    for (; srcItr != providedSockets.end(); ++srcItr)
    {
      if ((*meshData)->match((*srcItr)->getNamespace()))
      {
        CFLog(VERBOSE,"Allocating socket: " << (*srcItr)->getDataSocketName()
                      << " Storage: " << (*srcItr)->getDataSocketStorage()
                      << " Type: " << (*srcItr)->getDataSocketType()
                      << " Strategy: " << getName()
                      << "\n");
	
        (*srcItr)->allocate((*meshData)->getDataStorage(), 
			    (*srcItr)->getNamespace());
      }
    }
  }
  
  // check if any socket was not allocated
  VecSock::iterator srcItr = providedSockets.begin();
  for (; srcItr != providedSockets.end(); ++srcItr)
  {
    CF_DEBUG_OBJ ( (*srcItr)->getName() );
    if ( !(*srcItr)->isAllocated() )
    {
      throw BadMatchException (FromHere(),"Socket : " + (*srcItr)->getName()
            + " has invalid Namespace : " + (*srcItr)->getNamespace()
            + " in strategy: " + getName()
            + "\n" );
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NumericalStrategy::deallocateStrategySockets()
{
  vector<Common::SafePtr<BaseDataSocketSource> >::iterator srcItr;

  vector<Common::SafePtr<BaseDataSocketSource> > providedSockets = providesSockets();
  for (srcItr = providedSockets.begin(); srcItr != providedSockets.end(); ++srcItr) {

  CFLog(VERBOSE,"Deallocating socket: " << (*srcItr)->getDataSocketName()
                << " StorageType: " << (*srcItr)->getDataSocketStorage()
                << " Type: " << (*srcItr)->getDataSocketType()
                << "\n");

      (*srcItr)->deallocate();
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > NumericalStrategy::providesSockets()
{
  return std::vector<Common::SafePtr<BaseDataSocketSource> >();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > NumericalStrategy::needsSockets()
{
  return std::vector<Common::SafePtr<BaseDataSocketSink> >();
}

//////////////////////////////////////////////////////////////////////////////

void NumericalStrategy::setup()
{
  SetupObject::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NumericalStrategy::unsetup()
{
  SetupObject::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NumericalStrategy::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure (args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
