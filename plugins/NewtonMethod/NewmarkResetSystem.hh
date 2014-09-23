// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewmarkResetSystem_hh
#define COOLFluiD_Numerics_NewtonMethod_NewmarkResetSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "ResetSystem.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand of the NewtonIterator method
  /// used for resetting the linear system matrix, and right hand side vector
  /// to zero.
class NewmarkResetSystem : public ResetSystem {
public:

  /// Constructor.
  explicit NewmarkResetSystem(std::string name) :
    ResetSystem(name),
    socket_states("states"),
    socket_pastStates("pastStates")
  {
  }

  /// Destructor.
  ~NewmarkResetSystem()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  // socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // socket for pastStates
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

}; // class NewmarkResetSystem

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_NewmarkResetSystem_hh

