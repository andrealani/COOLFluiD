// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_TwoLayerUpdateSol_hh
#define COOLFluiD_Numerics_ForwardEuler_TwoLayerUpdateSol_hh

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
  /// @author Thomas Wuilbaut
class  ForwardEuler_API TwoLayerUpdateSol : public FwdEulerCom {
public:

  /// Constructor.
  explicit TwoLayerUpdateSol(const std::string& name) :
    FwdEulerCom(name),
    socket_pastStates("pastStates"),
    socket_interStates("interStates"),
    socket_rhs("rhs"),
    socket_interRhs("interRhs"),
    socket_updateCoeff("updateCoeff"),
    socket_interUpdateCoeff("interUpdateCoeff"),
    socket_states("states")
  {
  }

  /// Destructor.
  ~TwoLayerUpdateSol()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for PastStates's
  Framework::DataSocketSink<
                              Framework::State*> socket_pastStates;

  /// socket for InterStates's
  Framework::DataSocketSink<
                              Framework::State*> socket_interStates;

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

  /// socket for interRhs
  Framework::DataSocketSink<
                              CFreal> socket_interRhs;

  /// socket for updateCoeff
  Framework::DataSocketSink<
                              CFreal> socket_updateCoeff;

  /// socket for interUpdateCoeff
  Framework::DataSocketSink<
                              CFreal> socket_interUpdateCoeff;

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class TwoLayerUpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_TwoLayerUpdateSol_hh

