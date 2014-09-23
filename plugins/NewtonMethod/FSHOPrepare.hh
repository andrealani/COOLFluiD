// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_FSHOPrepare_hh
#define COOLFluiD_Numerics_NewtonMethod_FSHOPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
  /// @author Thomas Wuilbaut

class FSHOPrepare : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit FSHOPrepare(std::string name) :
    NewtonIteratorCom(name),
    socket_states("states"),
    socket_pastStates("pastStates"),
    socket_interStates("interStates")
  {
  }

  /// Destructor.
  ~FSHOPrepare()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // handle to rhs
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  // handle to update coefficient
  Framework::DataSocketSink<Framework::State*> socket_interStates;

}; // class FSHOPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_FSHOPrepare_hh

