// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/BadValueException.hh"

#include "Framework/Method.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/ConsistencyException.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MethodData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void Method::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Namespace","Namespace of this Methods.");
   options.addConfigOption< std::vector<std::string> >("CommandGroups","Names of the CommandGroups to be created.");
}

//////////////////////////////////////////////////////////////////////////////

Method::Method(const std::string& name)
  : Common::OwnedObject(),
    ConfigObject(name),
    SetupObject(),
    Common::NullableObject(),
    NamespaceMember(),
    m_commands(),
    m_groupNames(),
    m_groups(),
    m_isNonRootMethod(false)
{
  addConfigOptionsTo(this);

  m_groupNames = std::vector<std::string>();
  setParameter("CommandGroups",&m_groupNames);

  setParameter("Namespace",&m_namespace);
}

//////////////////////////////////////////////////////////////////////////////

Method::~Method()
{
  if (isSetup()) {
    unsetMethod();
  }
  
  for(CFuint i = 0; i < m_groups.size(); ++i) {
    deletePtr(m_groups[i]);
  }
}


//////////////////////////////////////////////////////////////////////////////

void Method::registActionListeners()
{
  CFAUTOTRACE;

#if defined CF_HAVE_BOOST_1_76 || defined CF_HAVE_BOOST_1_79
 create_signal ( "CF_ON_MESH_UPDATE" , "Resetup method when mesh updates" )->connect( boost::bind ( &Method::resetup, this, std::placeholders::_1 ) );
#else
  create_signal ( "CF_ON_MESH_UPDATE" , "Resetup method when mesh updates" )->connect( boost::bind ( &Method::resetup, this, _1 ) );
 #endif 
}

//////////////////////////////////////////////////////////////////////////////

void Method::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "Method::configure() for [" << getName() << "] => start\n");
  
  ConfigObject::configure(args);
  
  // configuring the socket namespace to which it belongs
  // this has precedence over the parent socket namespace
  if( m_namespace != NamespaceMember::defaultNamespace() )
  {
    setSelfNamespace(m_namespace);
  }
  else
  {
    std::vector<Common::SafePtr<Namespace> > allNsp = 
      NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getAllNamespaces();
    std::string defaultNamespace = allNsp[0]->getName();

    setParentNamespace(defaultNamespace);
  }
  
  CFLog(VERBOSE, "Configuring Method [" << getName() << "] in the Namespace: [" << getNamespace() << "]\n");
  
  SafePtr<MethodData> data = getMethodData();
  if (data.isNotNull())
  {
    data->setParentNamespace(getNamespace());
  }

  configureCommandGroups ( args );
  
  CFLog(VERBOSE, "Method::configure() for [" << getName() << "] => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void Method::configureCommandGroups ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Method::configureCommandGroups() for [" << getName() << "] => start\n");
  
  // configure the CommandGroups of this Method
  m_groups.resize(m_groupNames.size());

  CFLog(VERBOSE,"Method: " << getName() << " Number of CommandGroups: " << m_groups.size() << "\n");

  for (CFuint i = 0; i < m_groups.size(); ++i) {
    m_groups[i] = new CommandGroup(m_groupNames[i]);
    configureNested(m_groups[i], args);
  }
  
  CFLog(VERBOSE, "Method::configureCommandGroups() for [" << getName() << "] => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void Method::setMethod()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "Method::setMethod() for [" << getName() << "] => start\n");
  
  cf_assert(isConfigured());
  cf_assert(!isSetup());

  pushNamespace();
  

  // setup parent class
  SetupObject::setup();


  // setup method shared data
  Common::SafePtr<Framework::MethodData> mdata = getMethodData();
  if ( mdata.isNotNull() )
  {
    mdata->setup();
  }

  // setup derived classes
  this->setMethodImpl();

  popNamespace();
  
  CFLog(VERBOSE, "Method::setMethod() for [" << getName() << "] => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void Method::unsetMethod()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Method::unsetMethod() for [" << getName() << "] => start\n");
  
  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();
  
  // unsetup derived classes
  this->unsetMethodImpl();
  
  // unsetup method shared data
  Common::SafePtr<Framework::MethodData> mdata = getMethodData();
  if ( mdata.isNotNull() ) {
    mdata->unsetup();
  }
  
  // unsetup parent class
  SetupObject::unsetup();
  
  popNamespace(); 
  
  CFLog(VERBOSE, "Method::unsetMethod() for [" << getName() << "] => end\n");
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t
Method::resetup ( Common::Signal::arg_t input )
{
  CFAUTOTRACE;
  
  Method::unsetMethod();
  Method::setMethod();

  return Common::Signal::return_t();

}

//////////////////////////////////////////////////////////////////////////////

void Method::setupCommandsAndStrategies()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Method::setupCommandsAndStrategies() for [" << getName() << "] => start\n");
  
  vector< SafePtr<NumericalCommand> > commands = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comItr;
  for (comItr = commands.begin(); comItr != commands.end(); ++comItr)
  {
    SafePtr<NumericalCommand>  com = (*comItr);
    cf_assert(com != CFNULL);
    
    CFLog(VERBOSE,"NumericalCommand name : " << (*comItr)->getName() << " START\n");
    if (!com->isSetup()) com->setup();
    CFLog(VERBOSE,"NumericalCommand name : " << (*comItr)->getName() << " END\n");
  }
  
  vector<Common::SafePtr<Framework::NumericalStrategy> > strategies = getStrategyList();
  vector<Common::SafePtr<Framework::NumericalStrategy> >::iterator strtItr;
  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr)
  {
    Common::SafePtr<Framework::NumericalStrategy> st = (*strtItr);
    cf_assert(st.isNotNull());
    
    CFLog(VERBOSE,"Strategy name : " << (*strtItr)->getName() << " START\n");
    if (!st->isSetup()) st->setup();
    CFLog(VERBOSE,"Strategy name : " << (*strtItr)->getName() << " END\n");
  }
  
  CFLog(VERBOSE, "Method::setupCommandsAndStrategies() for [" << getName() << "] => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void Method::unsetupCommandsAndStrategies()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "Method::unsetupCommandsAndStrategies() for [" << getName() << "] => start\n");
  
  vector< SafePtr<NumericalCommand> > commands = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comItr;
  for (comItr = commands.begin(); comItr != commands.end(); ++comItr)
  {
     SafePtr<NumericalCommand>  com = (*comItr);
    cf_assert(com != CFNULL);

    CFLog(DEBUG_MIN,"NumericalCommand name : " << (*comItr)->getName() << "\n");
    if (com->isSetup()) com->unsetup();
  }

  vector<Common::SafePtr<Framework::NumericalStrategy> > strategies = getStrategyList();
  vector<Common::SafePtr<Framework::NumericalStrategy> >::iterator strtItr;
  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr)
  {
    Common::SafePtr<Framework::NumericalStrategy> st = (*strtItr);
    cf_assert(st.isNotNull());
    
    CFLog(VERBOSE, "Method::unsetupCommandsAndStrategies() for Strategy [" << (*strtItr)->getPolymorphicTypeName() << ", " << (*strtItr)->getName() << "] \n");
    if (st->isSetup()) {
      st->unsetup();
    }
    else {
      CFLog(WARN, "Method::unsetupCommandsAndStrategies() => Strategy [" << (*strtItr)->getPolymorphicTypeName() << ", " << (*strtItr)->getName() 
	    << "] is NOT properly setup()\n");
    }
  }
  
  CFLog(VERBOSE, "Method::unsetupCommandsAndStrategies() for [" << getName() << "] => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void Method::setCommandGroups()
{
   CFAUTOTRACE;

   CFLog(VERBOSE, "Method::setCommandGroups() for [" << getName() << "] => start\n");

   cf_assert(isConfigured());

  // First assign the CommandGroup Name corresponding to each NumericalCommand
  //  CFout << "Assigning CommandGroup Names" << "\n";
  vector<CommandGroup*>::iterator groupItr;
  for (groupItr = m_groups.begin(); groupItr != m_groups.end(); ++groupItr) {
    cf_assert((*groupItr) != CFNULL);
    CFLogDebugMin( "CommandGroup name = " << (*groupItr)->getName() << "\n");
    //CFout << "CommandGroup name = " << (*groupItr)->getName() << "\n";
    const vector<std::string> comNames = (*groupItr)->getComNames();

    vector< SafePtr<NumericalCommand> > commands = getCommandList();
    vector< SafePtr<NumericalCommand> >::iterator comItr;
    for (comItr = commands.begin(); comItr != commands.end(); ++comItr) {
      cf_assert((*comItr) != CFNULL);

      CFLog(DEBUG_MIN,"NumericalCommand name = " << (*comItr)->getName() << "\n");
      //CFout << "     NumericalCommand name = " << (*comItr)->getName() << "\n";
      const std::string commandName = (*comItr)->getName();

      // Loop over all Groups to check to which CommandGroup they belong
      for(CFuint i=0;i < comNames.size();++i)
        {
        if(comNames[i] == commandName){
          (*comItr)->setCommandGroupName((*groupItr)->getName());
          CFLog(DEBUG_MIN,"Setting the CommandGroup name to " << (*groupItr)->getName() << "\n");
        //  CFout << "     Setting the CommandGroup name to " << (*groupItr)->getName() << "\n";
          }
        }
    }
  }

  //CFout << "Assigning CommandGroup to each numerical command" << "\n";

  // Then assign the CommandGroup to each NumericalCommand
  vector< SafePtr<NumericalCommand> > commands = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comItr;
  for (comItr = commands.begin(); comItr != commands.end(); ++comItr) {
    cf_assert((*comItr) != CFNULL);

    CFLog(DEBUG_MIN,"NumericalCommand name = " << (*comItr)->getName() << "\n");
    //CFout << "     NumericalCommand name = " << (*comItr)->getName() << "\n";
    const std::string comGroupName = (*comItr)->getCommandGroupName();
    // Loop over all Groups to check to which CommandGroup they belong
    for(CFuint i=0;i < m_groupNames.size();++i)
    {
      if(m_groupNames[i] == comGroupName){
      //  CFout << "     Group Name: " << comGroupName << "\n";
        (*comItr)->setCommandGroup(m_groups[i]);
        }
    }
  } 
  
  CFLog(VERBOSE, "Method::setCommandGroups() for [" << getName() << "] => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void Method::setSocketNamespaces()
{
  CFAUTOTRACE;

  vector<Common::SafePtr<DataSocket> > allSockets = this->getAllSockets();

  vector<Common::SafePtr<DataSocket> >::iterator sck = allSockets.begin();
  for(;sck!=allSockets.end();++sck)
  {
    cf_assert(sck->isNotNull());
    (*sck)->setParentNamespace(getNamespace());
  }
}

//////////////////////////////////////////////////////////////////////////////

void Method::allocateMethodSockets()
{
  CFAUTOTRACE;
#if 0
  vector< SafePtr<NumericalCommand> > comList = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comitr;
  for (comitr = comList.begin(); comitr != comList.end(); ++comitr)
    (*comitr)->allocateCommandSockets();

  vector<Common::SafePtr<NumericalStrategy> > strategies = getStrategyList();
  vector<Common::SafePtr<NumericalStrategy> >::iterator strtItr;
  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr)
    (*strtItr)->allocateStrategySockets();
#endif
}

//////////////////////////////////////////////////////////////////////////////

void Method::deallocateMethodSockets()
{
  CFAUTOTRACE;
#if 0

  vector< SafePtr<NumericalCommand> > comList = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comitr;
  for (comitr = comList.begin(); comitr != comList.end(); ++comitr)
    (*comitr)->deallocateCommandSockets();

  vector<Common::SafePtr<NumericalStrategy> > strategies = getStrategyList();
  vector<Common::SafePtr<NumericalStrategy> >::iterator strtItr;
  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr)
    (*strtItr)->deallocateStrategySockets();

#endif
}

//////////////////////////////////////////////////////////////////////////////

void Method::checkAllSockets()
{
  CFAUTOTRACE;

  vector<Common::SafePtr<BaseDataSocketSink> > neededSockets =
    this->getNeededSockets();

  vector<Common::SafePtr<BaseDataSocketSink> >::iterator sink = neededSockets.begin();
  for(;sink!=neededSockets.end();++sink)
  {
    CFLog(DEBUG_MIN, " Method::checkAllSockets() => checking DataSocket connection in method " << getName()
            << " Name: "    << (*sink)->getDataSocketName()
            << " Storage: " << (*sink)->getDataSocketStorage()
            << " Type: "    << (*sink)->getDataSocketType()
            << " Namespace: " << (*sink)->getNamespace() << "\n\n");

    if((*sink)->isEssential() && !(*sink)->isConnected())
    {
      ostringstream msg;
      msg << " Essential DataSocket was not connected in method " + getName()
          << " Name: "    << (*sink)->getDataSocketName()
          << " Storage: " << (*sink)->getDataSocketStorage()
          << " Type: "    << (*sink)->getDataSocketType()
          << " Namespace: " << (*sink)->getNamespace();
      throw ConsistencyException (FromHere(),msg.str());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<DataSocket> >  Method::getAllSockets()
{
  CFAUTOTRACE;

  std::vector<Common::SafePtr<DataSocket> > result;

  vector< SafePtr<NumericalCommand> > commands = getCommandList();
  vector<Common::SafePtr<NumericalStrategy> > strategies = getStrategyList();

  vector< SafePtr<NumericalCommand> >::iterator comItr;
  vector< SafePtr<NumericalStrategy> >::iterator strtItr;
  vector< SafePtr<BaseDataSocketSink> >::iterator sinkItr;
  vector< SafePtr<BaseDataSocketSource> >::iterator srcItr;

  for (comItr = commands.begin(); comItr != commands.end(); ++comItr) {

    // put the Sink Sockets
    std::vector<Common::SafePtr<BaseDataSocketSink> > needSockets = (*comItr)->needsSockets();
    for (sinkItr = needSockets.begin(); sinkItr != needSockets.end(); ++sinkItr) {
      result.push_back(sinkItr->d_castTo<DataSocket>());
    }

    // put the Source Sockets
    std::vector<Common::SafePtr<BaseDataSocketSource> > providedSockets = (*comItr)->providesSockets();
    for (srcItr = providedSockets.begin(); srcItr != providedSockets.end(); ++srcItr) {
      result.push_back(srcItr->d_castTo<DataSocket>());
    }
  }

  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr) {

    // put the Sink Sockets
    std::vector<Common::SafePtr<BaseDataSocketSink> > needSockets = (*strtItr)->needsSockets();
    for (sinkItr = needSockets.begin(); sinkItr != needSockets.end(); ++sinkItr) {
      result.push_back(sinkItr->d_castTo<DataSocket>());
    }

    // put the Source Sockets
    std::vector<Common::SafePtr<BaseDataSocketSource> > providedSockets = (*strtItr)->providesSockets();
    for (srcItr = providedSockets.begin(); srcItr != providedSockets.end(); ++srcItr) {
      result.push_back(srcItr->d_castTo<DataSocket>());
    }
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >  Method::getNeededSockets()
{
  CFAUTOTRACE;

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  CFLog(DEBUG_MIN, "Method::getNeededSockets() => " << getName() << " Commands: ["); 
  vector< SafePtr<NumericalCommand> > commands = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comItr;
  for (comItr = commands.begin(); comItr != commands.end(); ++comItr) {
    CFLog(DEBUG_MIN, (*comItr)->getName() << " ");
    std::vector<Common::SafePtr<BaseDataSocketSink> > comSockets = (*comItr)->needsSockets();
    copy(comSockets.begin(),comSockets.end(),back_inserter(result));
  }
  CFLog(DEBUG_MIN,"]\n");

   CFLog(DEBUG_MIN, "Method::getNeededSockets() => " << getName() << " Strategies: [");
  vector<Common::SafePtr<NumericalStrategy> > strategies = getStrategyList();
  vector<Common::SafePtr<NumericalStrategy> >::iterator strtItr;
  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr) {
    CFLog(DEBUG_MIN, (*strtItr)->getName() << " ");
    std::vector<Common::SafePtr<BaseDataSocketSink> > stratSockets = (*strtItr)->needsSockets();
    copy(stratSockets.begin(),stratSockets.end(),back_inserter(result));
  }
  CFLog(DEBUG_MIN,"]\n");

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > Method::getProvidedSockets()
{
  CFTRACEBEGIN;

  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  // get sockets from all commands
  vector< SafePtr<NumericalCommand> > commands = getCommandList();
  vector< SafePtr<NumericalCommand> >::iterator comItr;
  for (comItr = commands.begin(); comItr != commands.end(); ++comItr)
  {
    std::vector<Common::SafePtr<BaseDataSocketSource> > comSockets = (*comItr)->providesSockets();
    copy(comSockets.begin(),comSockets.end(),back_inserter(result));
  }

  // get sockets from all stategies
  vector<Common::SafePtr<NumericalStrategy> > strategies = getStrategyList();
  vector<Common::SafePtr<NumericalStrategy> >::iterator strtItr;
  for (strtItr = strategies.begin(); strtItr != strategies.end(); ++strtItr)
  {
    std::vector<Common::SafePtr<BaseDataSocketSource> > stratSockets = (*strtItr)->providesSockets();
    copy(stratSockets.begin(),stratSockets.end(),back_inserter(result));
  }

  CFTRACEEND;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<CommandGroup> > Method::getCommandGroups()
{
  std::vector<Common::SafePtr<CommandGroup> > result;

  std::vector<CommandGroup*>::iterator itr = m_groups.begin();
  for(; itr != m_groups.end(); ++itr) {
    result.push_back(*itr);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Method::pushNamespace()
{
  CFAUTOTRACE;
  
  NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).
    pushNamespace(getNamespace());
}

//////////////////////////////////////////////////////////////////////////////

void Method::popNamespace()
{
  CFAUTOTRACE;

  SafePtr<Namespace> ptr = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).popNamespace();
  cf_assert(ptr->getName() == getNamespace());
}

//////////////////////////////////////////////////////////////////////////////

QualifiedName Method::QName() const
{
  return QualifiedName (getNamespace(), getName());
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr<NumericalCommand> >
Method::getCommandList() const
{
  return m_commands;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::NumericalStrategy> >
Method::getStrategyList() const
{
  Common::SafePtr<Framework::MethodData> mdata = getMethodData();
  if ( mdata.isNotNull() ) // null methods might not have data
    return mdata->getStrategies();
  else // empty vector
    return std::vector<Common::SafePtr<Framework::NumericalStrategy> >();
}

//////////////////////////////////////////////////////////////////////////////

void Method::executeCommands(const vector<std::string>& comNames)
{
  pushNamespace();

  vector< Common::SafePtr<NumericalCommand> > comList = getCommandList();
  for (CFuint i = 0; i < comNames.size(); ++i)
  {
    bool nameFound = false;
    for (CFuint j = 0; j < comList.size(); ++j)
    {
      if (comNames[i] == comList[j]->getName())
      {
        comList[j]->execute();
        nameFound = true;
        break;
      }
    }
    if (!nameFound)
    {
      ostringstream msg;
      msg << "Calling non-existing command [" << comNames[i] << "]\n" << std::endl;
      throw BadValueException ( FromHere(), msg.str() );
    }
  }

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
