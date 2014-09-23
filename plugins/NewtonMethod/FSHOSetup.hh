// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_FSHOSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_FSHOSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action
/// to be executed in order to setup the MeshData.
/// @author Thomas Wuilbaut
class FSHOSetup : public StdSetup {
public:

  /// Constructor.
  explicit FSHOSetup(std::string name);

  /// Destructor.
  ~FSHOSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Execute Processing actions
  void execute();

private:

  /// socket for InterStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_interStates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_FSHOSetup_hh
