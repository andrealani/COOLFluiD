// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptySpaceMethod_StdSetup_hh
#define COOLFluiD_EmptySpaceMethod_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "EmptySpaceMethod/EmptySolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to setup the empty method
/// @author Tiago Quintino
/// @author Pedro Maciel
class StdSetup : public EmptySolverCom {

public: // functions

  /// Constructor
  explicit StdSetup(const std::string& name);

  /// Destructor
  ~StdSetup();

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

  /// Node to state ID map
  std::vector< CFuint > m_nodeIDToStateID;

protected: // data

  /// socket with proxy to be able to use the data handle of nodal states
  /// uniformly independently from the actual storage type being RealVector or
  /// State*
  Framework::DataSocketSource< Framework::ProxyDofIterator< RealVector >* >
    socket_nstatesProxy;

  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

};  // class StdSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptySpaceMethod_StdSetup_hh

