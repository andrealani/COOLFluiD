#ifndef COOLFluiD_network_NetworkFrameType_h
#define COOLFluiD_network_NetworkFrameType_h

/////////////////////////////////////////////////////////////////////////////

#include "ClientServer/network/NetworkAPI.hh"

namespace COOLFluiD
{
  namespace network
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Defines all frame types allowed in the network m_protocol.
    
    enum NetworkFrameType 
    {
      /// @brief Type id used to indicate that a frame has no type (i.e. it 
      /// does not respect the m_protocol).
      NETWORK_NO_TYPE,
      
      /// @brief Type id for the frame root tag.
      NETWORK_FRAME_ROOT,
      
      /// @brief Type id for "Error" frame.
      NETWORK_ERROR,
      
      /// @brief Type id for "Get tree" frame.
      NETWORK_GET_TREE,
      
      /// @brief Type id for "Message" frame.
      NETWORK_MESSAGE,
      
      /// @brief Type id for "Modify node" frame.
      NETWORK_MODIFY_NODE,
      
      /// @brief Type id for "Tree" frame.
      NETWORK_TREE,
      
      /// @brief Type id for "Add node" frame.
      NETWORK_ADD_NODE,
      
      /// @brief Type id for "Delete node" frame.
      NETWORK_DELETE_NODE,
      
      /// @brief Type id for "Rename" frame.
      NETWORK_RENAME_NODE,
      
      /// @brief Type id for "Get abstract types" frame.
      NETWORK_GET_ABSTRACT_TYPES,
      
      /// @brief Type id for "Get concrete types" frame.
      NETWORK_GET_CONCRETE_TYPES,
      
      /// @brief Type id for "Abstract types" frame.
      NETWORK_ABSTRACT_TYPES,
      
      /// @brief Type id for "Concrete types" frame.
      NETWORK_CONCRETE_TYPES,
      
      /// @brief Type id for "Get files list" frame.
      NETWORK_GET_FILES_LIST,
      
      /// @brief Type id for "Files list" frame.
      NETWORK_FILE_LIST,
      
      /// @brief Type id for "Open file" frame.
      NETWORK_OPEN_FILE,
      
      /// @brief Type id for "Run simulation" frame.
      NETWORK_RUN_SIMULATION,
      
      /// @brief Type id for "ACK" frame.
      NETWORK_ACK,
      
      /// @brief Type id for "NACK" frame.
      NETWORK_NACK,
      
      /// @brief Type id for "Shutdown server" frame.
      NETWORK_SHUTDOWN_SERVER,
      
      /// @brief Type id for "Simulation running" frame.
      NETWORK_SIMULATION_RUNNING,
      
      /// @brief Type id for "Open directory" frame.
      NETWORK_READ_DIR,
      
      /// @brief Type id for "Directory contents" frame.
      NETWORK_DIR_CONTENTS,
      
      /// @brief Type id for "Close file" frame.
      NETWORK_CLOSE_FILE,
      
      /// @brief Type id for "Configure" frame.
      NETWORK_CONFIGURE,
      
      /// @brief Type id for "Create directory" frame.
      NETWORK_CREATE_DIR,
      
      /// @brief Type id for "Save configuration" frame.
      NETWORK_SAVE_CONFIG,
      
      /// @brief Type id for "Get host list" frame.
      NETWORK_GET_HOST_LIST,
      
      /// @brief Type id for "Host list" frame.
      NETWORK_HOST_LIST,
      
      /// @brief Type id for "Stop configuration" frame.
      NETWORK_STOP_SIMULATION,
      
      /// @brief Type id for "Activate simulation" frame.
      NETWORK_ACTIVATE_SIMULATION,
      
      /// @brief Type id for "Deactivate simulation" frame.
      NETWORK_DEACTIVATE_SIMULATION,
      
      /// @brief Type id for "Simulation status" frame.
      NETWORK_SIMULATION_STATUS,
      
      NETWORK_SUBSYSTEM_LIST,
      
      NETWORK_GET_SUBSYSTEM_LIST
      
    }; // enum FrameTypes
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace network
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_NetworkFrameType_h
