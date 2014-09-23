// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_StdPrepare_hh
#define COOLFluiD_Numerics_ForwardEuler_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "FwdEulerData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class  ForwardEuler_API StdPrepare : public FwdEulerCom {
public:

  /// Constructor.
  explicit StdPrepare(std::string name) :
    FwdEulerCom(name),
    socket_states("states",false),
    socket_pastStates("pastStates",false)
  {
  }

  /// Destructor.
  ~StdPrepare()
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

  // handle to past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

}; // class StdPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_StdPrepare_hh
