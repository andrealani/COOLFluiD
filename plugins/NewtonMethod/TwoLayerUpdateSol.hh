// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_TwoLayerUpdateSol_hh
#define COOLFluiD_Numerics_NewtonMethod_TwoLayerUpdateSol_hh

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

class TwoLayerUpdateSol : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit TwoLayerUpdateSol(const std::string& name) :
    NewtonIteratorCom(name),
    socket_states("states"),
    socket_interStates("interStates"),
    socket_rhs("rhs"),
    socket_interRhs("interRhs"),
    socket_updateCoeff("updateCoeff"),
    socket_interUpdateCoeff("interUpdateCoeff")
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

protected:

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to intermediate states
  Framework::DataSocketSink< Framework::State*> socket_interStates;

  /// handle to rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

  /// handle to intermediate rhs
  Framework::DataSocketSink<CFreal> socket_interRhs;

  // handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  // handle to intermediate update coefficient
  Framework::DataSocketSink<CFreal> socket_interUpdateCoeff;

}; // class TwoLayerUpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_TwoLayerUpdateSol_hh
