// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Config/BadMatchException.hh"
#include "Common/CFLog.hh"
#include "Framework/NumericalCommand.hh"
#include "Framework/BaseDataSocketSource.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/MeshData.hh"
#include "Framework/ConsistencyException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("ComGroup","Name of the CommandGroup to which this Command belongs.");
  options.addConfigOption< std::vector<std::string> >("applyTRS","Names of the TRS's to which this Command applies.");
  options.addConfigOption< CFuint >
    ("ProcessRate", "Rate at which the processing has to be executed.");
}

//////////////////////////////////////////////////////////////////////////////

NumericalCommand::NumericalCommand(const std::string& name)
  : Common::OwnedObject(),
    ConfigObject(name),
    m_iTrs(0),
    m_trsList(),
    m_group(CFNULL)
{
  addConfigOptionsTo(this);
  m_trsNames = std::vector<std::string>();
  setParameter("applyTRS",&m_trsNames);

  m_groupName = "";
  setParameter("ComGroup",&m_groupName);

  m_processRate = 1000000000;
  setParameter("ProcessRate", &m_processRate);
}

//////////////////////////////////////////////////////////////////////////////

NumericalCommand::~NumericalCommand()
{
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::execute()
{
  CFuint nbTrs = m_trsList.size();
  CFLogDebugMed("Command: " << getName() << " will be executed in " << nbTrs << " TRSs" << "\n");
  for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs) {
    CFLogDebugMed("Command: " << getName() << " applying on TRS: " << (m_trsList[iTrs])->getName() << "\n");
    setCurrentTrsID(iTrs);
    executeOnTrs();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::setCommandGroup(Common::SafePtr<CommandGroup> commandGroup)
{
  m_group = commandGroup;
  std::vector<std::string> trsNames = m_group->getTrsNames();
  m_trsNames.resize(trsNames.size());
  m_trsNames = trsNames;
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::configureNestedSockets( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure the source sockets
  vector<Common::SafePtr<BaseDataSocketSource> >::iterator srcItr;
  vector<Common::SafePtr<BaseDataSocketSource> > providedSockets = providesSockets();
  for (srcItr = providedSockets.begin(); srcItr != providedSockets.end(); ++srcItr) {
      cf_assert(srcItr->isNotNull());
      configureNested(*(*srcItr), args);
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

void NumericalCommand::allocateCommandSockets()
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
    VecSock::iterator srcItr = providedSockets.begin();
    for (; srcItr != providedSockets.end(); ++srcItr)
    {
      if ((*meshData)->match((*srcItr)->getNamespace()))
      {
        CFLog(VERBOSE,"Allocating socket: " << (*srcItr)->getDataSocketName()
                      << " Storage: " << (*srcItr)->getDataSocketStorage()
                      << " Type: " << (*srcItr)->getDataSocketType()
                      << " Command: " << getName()
                      << "\n");

        (*srcItr)->allocate((*meshData)->getDataStorage(), (*srcItr)->getNamespace());
      }
    }
  }

  // check if any socket was not allocated
  VecSock::iterator srcItr = providedSockets.begin();
  for (; srcItr != providedSockets.end(); ++srcItr)
  {
    if (!((*srcItr)->isAllocated()))
    {
      throw BadMatchException (FromHere(),"Socket : " + (*srcItr)->getName()
            + " in Namespace : " + (*srcItr)->getNamespace()
            + " in command: " + getName()
            + " was not allocated.\n" );
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::deallocateCommandSockets()
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

std::vector<Common::SafePtr<BaseDataSocketSource> > NumericalCommand::providesSockets()
{
  return std::vector<Common::SafePtr<BaseDataSocketSource> >();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > NumericalCommand::needsSockets()
{
  return std::vector<Common::SafePtr<BaseDataSocketSink> >();
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::setup()
{
  SetupObject::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::unsetup()
{
  SetupObject::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::executeOnTrs()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<TopologicalRegionSet> NumericalCommand::getCurrentTRS() const
{
  if (m_trsList.empty())
    throw ConsistencyException (FromHere(),"Calling getCurrentTRS() on a command which has no TRS set");

  return m_trsList[m_iTrs];
}

//////////////////////////////////////////////////////////////////////////////

void NumericalCommand::setCurrentTrsID(const CFuint iTrs)
{
  m_iTrs = iTrs;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr<TopologicalRegionSet> >&
NumericalCommand::getTrsList()
{
  return m_trsList;
}

//////////////////////////////////////////////////////////////////////////////

const std::string NumericalCommand::getTrsName(const CFuint iTrs)
{
  if (m_trsNames.empty())
    throw ConsistencyException (FromHere(),"Calling getTrsName() on a command which has no TRS assigned");

  return m_trsNames[iTrs];
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
