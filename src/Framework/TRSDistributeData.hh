// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TRSDistributeData_hh
#define COOLFluiD_Framework_TRSDistributeData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////
    
class Framework_API TRSDistributeData {
public:
  
  std::vector<CFreal> trsNodes;
  std::vector<CFint>  trsNodeConn;
  std::vector<CFint>  trsStateConn;
  std::vector<CFuint> trsNbNodesInFace;
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TRGeoConn_hh
