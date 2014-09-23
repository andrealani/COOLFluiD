// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Pardiso_StdSolveSys_hh
#define COOLFluiD_Pardiso_StdSolveSys_hh

#include "Pardiso/PardisoData.hh"
#include "Framework/DataSocketSink.hh"

namespace COOLFluiD {
  namespace Pardiso {

/// This is a standard command to solve the linear system using PARDISO
class StdSolveSys : public PardisoCom {

public:

  /// Constructor
  explicit StdSolveSys(const std::string& name) :
    PardisoCom(name),
    socket_rhs("rhs") {}

  /// Destructor
  virtual ~StdSolveSys() {}

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets() {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > result;
    result.push_back(&socket_rhs);
    return result;
  }


private: // data

  /// Handle to RHS
  Framework::DataSocketSink<CFreal > socket_rhs;

}; // class SolveSys


  } // namespace Pardiso
} // namespace COOLFluiD

#endif // COOLFluiD_Pardiso_StdSolveSys_hh

