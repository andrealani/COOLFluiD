#include <QString>

#include "ClientServer/network/NetworkFrameType.hh"
#include "ClientServer/network/NetworkProtocol.hh"

using namespace COOLFluiD::network;
using namespace COOLFluiD::Config;

NetworkProtocol::NetworkProtocol() 
: BuilderParserRules(NETWORK_NO_TYPE, NETWORK_FRAME_ROOT, "ClientServerXML")
{
  this->buildRules();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NetworkFrameType NetworkProtocol::convertToType(unsigned int type) const
{
  switch(type)
  {
    case NETWORK_FRAME_ROOT:
      return NETWORK_FRAME_ROOT;
      
    case NETWORK_ABSTRACT_TYPES:
      return NETWORK_ABSTRACT_TYPES;
      
    case NETWORK_ACK:
      return NETWORK_ACK;
      
    case NETWORK_ADD_NODE:
      return NETWORK_ADD_NODE;
      
    case NETWORK_CLOSE_FILE:
      return NETWORK_CLOSE_FILE;
      
    case NETWORK_CONCRETE_TYPES:
      return NETWORK_CONCRETE_TYPES;
      
    case NETWORK_CONFIGURE:
      return NETWORK_CONFIGURE;
      
    case NETWORK_CREATE_DIR:
      return NETWORK_CREATE_DIR;
      
    case NETWORK_DELETE_NODE:
      return NETWORK_DELETE_NODE;
      
    case NETWORK_DIR_CONTENTS:
      return NETWORK_DIR_CONTENTS;
      
    case NETWORK_ERROR:
      return NETWORK_ERROR;
      
    case NETWORK_GET_ABSTRACT_TYPES:
      return NETWORK_GET_ABSTRACT_TYPES;
      
    case NETWORK_GET_CONCRETE_TYPES:
      return NETWORK_GET_CONCRETE_TYPES;
      
    case NETWORK_GET_TREE:
      return NETWORK_GET_TREE;
      
    case NETWORK_MESSAGE:
      return NETWORK_MESSAGE;
      
    case NETWORK_MODIFY_NODE:
      return NETWORK_MODIFY_NODE;
      
    case NETWORK_NACK:
      return NETWORK_NACK;
      
    case NETWORK_OPEN_FILE:
      return NETWORK_OPEN_FILE;
      
    case NETWORK_READ_DIR:
      return NETWORK_READ_DIR;
      
    case NETWORK_RENAME_NODE:
      return NETWORK_RENAME_NODE;
      
    case NETWORK_RUN_SIMULATION:
      return NETWORK_RUN_SIMULATION;
      
    case NETWORK_SAVE_CONFIG:
      return NETWORK_SAVE_CONFIG;
      
    case NETWORK_SHUTDOWN_SERVER:
      return NETWORK_SHUTDOWN_SERVER;
      
    case NETWORK_SIMULATION_RUNNING:
      return NETWORK_SIMULATION_RUNNING;
      
    case NETWORK_TREE:
      return NETWORK_TREE;
      
    case NETWORK_STOP_SIMULATION:
      return NETWORK_STOP_SIMULATION;
      
    case NETWORK_ACTIVATE_SIMULATION:
      return NETWORK_ACTIVATE_SIMULATION;
      
    case NETWORK_DEACTIVATE_SIMULATION:
      return NETWORK_DEACTIVATE_SIMULATION;
      
    case NETWORK_SIMULATION_STATUS:
      return NETWORK_SIMULATION_STATUS;
      
    case NETWORK_SUBSYSTEM_LIST:
      return NETWORK_SUBSYSTEM_LIST;
      
    case NETWORK_GET_SUBSYSTEM_LIST:
      return NETWORK_GET_SUBSYSTEM_LIST;
      
    default:
      return NETWORK_NO_TYPE;
  } 
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NetworkFrameType NetworkProtocol::convertToType(const QString & typeName) const
{
  return this->convertToType(this->getFrameType(typeName.toStdString()));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void NetworkProtocol::buildRules()
{
  // Abstract types
  this->setTypeName(NETWORK_ABSTRACT_TYPES, "abstractTypes");
  this->addAttribute(NETWORK_ABSTRACT_TYPES, "typeName");
  this->addAttribute(NETWORK_ABSTRACT_TYPES, "typeList");
  
  // ACK
  this->setTypeName(NETWORK_ACK, "ack");
  this->addAttribute(NETWORK_ACK, "type");
  
  // Add node
  this->setTypeName(NETWORK_ADD_NODE, "addNode");
  this->addAttribute(NETWORK_ADD_NODE, "path");
  this->addAttribute(NETWORK_ADD_NODE, "type");
  this->addAttribute(NETWORK_ADD_NODE, "absType");
  
  // Close file
  this->setTypeName(NETWORK_CLOSE_FILE, "closeFile");
  
  // Concrete types
  this->setTypeName(NETWORK_CONCRETE_TYPES, "concreteTypes");
  this->addAttribute(NETWORK_CONCRETE_TYPES, "typeName");
  this->addAttribute(NETWORK_CONCRETE_TYPES, "typeList");
  
  // Configure
  this->setTypeName(NETWORK_CONFIGURE, "configure");
  this->setDataMandatoriness(NETWORK_CONFIGURE, MAND_MANDATORY);
  
  // Create directory
  this->setTypeName(NETWORK_CREATE_DIR, "createDirectory");
  this->addAttribute(NETWORK_CREATE_DIR, "path");
  this->addAttribute(NETWORK_CREATE_DIR, "dirName");
  
  // Delete node
  this->setTypeName(NETWORK_DELETE_NODE, "deleteNode");
  this->addAttribute(NETWORK_DELETE_NODE, "path");
  
  // Directory contents
  this->setTypeName(NETWORK_DIR_CONTENTS, "dirContents");
  this->addAttribute(NETWORK_DIR_CONTENTS, "dirPath");
  this->addAttribute(NETWORK_DIR_CONTENTS, "dirs");
  this->addAttribute(NETWORK_DIR_CONTENTS, "files", MAND_OPTIONAL);
  
  // Error
  this->setTypeName(NETWORK_ERROR, "error");
  this->addAttribute(NETWORK_ERROR, "value");
  
  // Get abstract types
  this->setTypeName(NETWORK_GET_ABSTRACT_TYPES, "getAbstractTypes");
  this->addAttribute(NETWORK_GET_ABSTRACT_TYPES, "typeName");
  
  // Get concrete types
  this->setTypeName(NETWORK_GET_CONCRETE_TYPES, "getConcreteTypes");
  this->addAttribute(NETWORK_GET_CONCRETE_TYPES, "typeName");
  
  // Get host list
  this->setTypeName(NETWORK_GET_HOST_LIST, "getHostList");
  
  // Get tree
  this->setTypeName(NETWORK_GET_TREE, "getTree");
  
  // Host list
  this->setTypeName(NETWORK_HOST_LIST, "hostList");
  this->setDataMandatoriness(NETWORK_HOST_LIST, MAND_OPTIONAL);
  
  // Message
  this->setTypeName(NETWORK_MESSAGE, "message");
  this->addAttribute(NETWORK_MESSAGE, "value");
  
  // Modify node
  this->setTypeName(NETWORK_MODIFY_NODE, "modifyNode");
  this->setDataMandatoriness(NETWORK_MODIFY_NODE, MAND_MANDATORY);
  
  // NACK
  this->setTypeName(NETWORK_NACK, "nack");
  this->addAttribute(NETWORK_NACK, "type");
  
  // Open file
  this->setTypeName(NETWORK_OPEN_FILE, "openFile");
  this->addAttribute(NETWORK_OPEN_FILE, "filename");
  
  // Read dir
  this->setTypeName(NETWORK_READ_DIR, "readDir");
  this->addAttribute(NETWORK_READ_DIR, "dirPath", MAND_OPTIONAL);
  this->addAttribute(NETWORK_READ_DIR, "includeFiles");
  this->addAttribute(NETWORK_READ_DIR, "extensions");
  this->addAttribute(NETWORK_READ_DIR, "includeNoExtension");
  
  // Rename node
  this->setTypeName(NETWORK_RENAME_NODE, "renameNode");
  this->addAttribute(NETWORK_RENAME_NODE, "path");
  this->addAttribute(NETWORK_RENAME_NODE, "newName");
  
  // Run simulation
  this->setTypeName(NETWORK_RUN_SIMULATION, "runSimulation");
  
  // Save configuration
  this->setTypeName(NETWORK_SAVE_CONFIG, "saveConfig");
  this->addAttribute(NETWORK_SAVE_CONFIG, "filename");
  this->setDataMandatoriness(NETWORK_SAVE_CONFIG, MAND_MANDATORY);
  
  // Shutdown server
  this->setTypeName(NETWORK_SHUTDOWN_SERVER, "shutdownServer");
  
  // Simulation running
  this->setTypeName(NETWORK_SIMULATION_RUNNING, "simulationRunning");
  
  // Tree
  this->setTypeName(NETWORK_TREE, "tree");
  this->setDataMandatoriness(NETWORK_TREE, MAND_MANDATORY);
  
  // Stop simulation
  this->setTypeName(NETWORK_STOP_SIMULATION, "stopSimulation");
  
  // Activate simulation
  this->setTypeName(NETWORK_ACTIVATE_SIMULATION, "activateSimulation");
  this->addAttribute(NETWORK_ACTIVATE_SIMULATION, "nbProcs");
  this->addAttribute(NETWORK_ACTIVATE_SIMULATION, "hosts", MAND_OPTIONAL);
  
  // Deactivate simulation
  this->setTypeName(NETWORK_DEACTIVATE_SIMULATION, "deactivateSimulation");
  
  // Simulation status
  this->setTypeName(NETWORK_SIMULATION_STATUS, "simulationStatus");
  this->addAttribute(NETWORK_SIMULATION_STATUS, "subSysName", MAND_OPTIONAL);
  this->addAttribute(NETWORK_SIMULATION_STATUS, "workerRank");
  this->addAttribute(NETWORK_SIMULATION_STATUS, "value");
  
  this->setTypeName(NETWORK_SUBSYSTEM_LIST, "subSystemList");
  this->addAttribute(NETWORK_SUBSYSTEM_LIST, "subSystems");
  
  this->setTypeName(NETWORK_GET_SUBSYSTEM_LIST, "getSubSystemList");
}
