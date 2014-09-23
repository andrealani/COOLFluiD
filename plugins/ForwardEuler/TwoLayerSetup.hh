// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_TwoLayerSetup_hh
#define COOLFluiD_Numerics_ForwardEuler_TwoLayerSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FwdEulerData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action
/// to be executed in order to setup the MeshData.
/// @author Thomas Wuilbaut
class  ForwardEuler_API TwoLayerSetup : public FwdEulerCom {
public:

  /// Constructor.
  explicit TwoLayerSetup(std::string name);

  /// Destructor.
  ~TwoLayerSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for PastStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_pastStates;

  /// socket for InterStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_interStates;

  /// socket for Rhs
  Framework::DataSocketSource<
                              CFreal> socket_rhs;

  /// socket for interRhs
  Framework::DataSocketSource<
                              CFreal> socket_interRhs;

  /// socket for updateCoeff
  Framework::DataSocketSource<
                              CFreal> socket_updateCoeff;

  /// socket for interUpdateCoeff
  Framework::DataSocketSource<
                              CFreal> socket_interUpdateCoeff;

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_TwoLayerSetup_hh

