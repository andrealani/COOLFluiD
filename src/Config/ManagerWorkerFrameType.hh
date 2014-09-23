// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_ManagerWorkerFrameType_hh
#define COOLFluiD_Config_ManagerWorkerFrameType_hh

#include "Config/Config.hh"

namespace COOLFluiD {
 
 namespace Config {
  
  enum ManagerWorkerFrameType
  {
   /// Type id used to indicate that a frame has no type (i.e. it 
   /// does not respect the protocol).
   MGR_WKR_NO_TYPE,
  
   /// Type id for the frame root tag.
   MGR_WKR_FRAME_ROOT,
  
   MGR_WKR_CONNECT,
   
   MGR_WKR_OPEN_PORT,
      
   MGR_WKR_PORT_NAME,
   
   MGR_WKR_CONFIGURE,
   
   MGR_WKR_SIMULATE,
   
   MGR_WKR_ACK,

   MGR_WKR_STRING,
   
   MGR_WKR_SET_SUBSYS,
   
   MGR_WKR_EXIT,
   
   MGR_WKR_STATUS,
   
   MGR_WKR_TREE,
   
   MGR_WKR_QUIT
     
  }; // enum ManagerWorkerFrameType
  
 } // namespace Config
} // namespace COOFluiD

#endif // COOLFluiD_Config_ManagerWorkerFrameType_hh