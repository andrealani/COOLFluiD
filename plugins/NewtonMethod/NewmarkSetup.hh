// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewmarkSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_NewmarkSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
  /// @author Thomas Wuilbaut
class NewmarkSetup : public StdSetup {
public:

  /// Constructor.
  explicit NewmarkSetup(std::string name);

  /// Destructor.
  ~NewmarkSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Execute Processing actions
  void execute();

protected:

  /// socket for PastStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_pastStatesD;

  /// socket for PastStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_pastStatesD2;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_NewmarkSetup_hh
