// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_LinearizedPrepare_hh
#define COOLFluiD_Numerics_NewtonMethod_LinearizedPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdPrepare.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class LinearizedPrepare : public StdPrepare {
public:

  /// Constructor.
  explicit LinearizedPrepare(std::string name) :
    StdPrepare(name),
    socket_linearizedStates("linearizedStates"),
    socket_linearizedGhostStates("linearizedGhostStates"),
    socket_pastGhostStates("pastGhostStates"),
    socket_gstates("gstates")
  {
  }

  /// Destructor.
  ~LinearizedPrepare()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  // handle to linearized states
  Framework::DataSocketSink<Framework::State*> socket_linearizedStates;

  // handle to linearized Ghost states
  Framework::DataSocketSink<Framework::State*> socket_linearizedGhostStates;

  // handle to past Ghost states
  Framework::DataSocketSink<Framework::State*> socket_pastGhostStates;

  // handle to Ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

}; // class LinearizedPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_LinearizedPrepare_hh

