// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_CrankNichLimUnSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_CrankNichLimUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CrankNichUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Thomas Wuilbaut
class CrankNichLimUnSetup : public CrankNichUnSetup {
public:

  /// Constructor.
  explicit CrankNichLimUnSetup(std::string name);

  /// Destructor.
  ~CrankNichLimUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected:

  /// socket for the past past states
  Framework::DataSocketSink<
                            Framework::State*> socket_pastPastStates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_CrankNichLimUnSetup_hh
