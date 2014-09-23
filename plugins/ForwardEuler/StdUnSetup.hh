// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_StdUnSetup_hh
#define COOLFluiD_Numerics_ForwardEuler_StdUnSetup_hh

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
/// @author Tiago Quintino
/// @author Andrea Lani
class  ForwardEuler_API StdUnSetup : public FwdEulerCom {
public:

  /// Constructor.
  explicit StdUnSetup(std::string name);

  /// Destructor.
  ~StdUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected:

  /// socket for PastStates's
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_StdUnSetup_hh

