// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_LinearizedUnSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_LinearizedUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Andrea Lani
class LinearizedUnSetup : public StdUnSetup {
public:

  /// Constructor.
  explicit LinearizedUnSetup(std::string name);

  /// Destructor.
  ~LinearizedUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected: //member data

  /// socket for the linearized states
  Framework::DataSocketSink<Framework::State*> socket_linearizedStates;


}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_LinearizedUnSetup_hh

