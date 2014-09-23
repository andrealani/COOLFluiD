// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_FSHOUnSetup_hh
#define COOLFluiD_Numerics_ForwardEuler_FSHOUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FwdEulerData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Nadege Villedieu
class  ForwardEuler_API FSHOUnSetup : public FwdEulerCom {
public:

  /// Constructor.
  explicit FSHOUnSetup(std::string name);

  /// Destructor.
  ~FSHOUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected:

  /// socket for PastStates's
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// socket for InterStates's
  Framework::DataSocketSink<Framework::State*> socket_interStates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_StdUnSetup_hh

