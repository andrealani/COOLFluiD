// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_TwoLayerSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_TwoLayerSetup_hh

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
class TwoLayerSetup : public StdSetup {
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

private:

  /// socket for InterStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_interStates;

  /// socket for interRhs
  Framework::DataSocketSource<
                              CFreal> socket_interRhs;

  /// socket for interUpdateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSource<
                              CFreal> socket_interUpdateCoeff;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_TwoLayerSetup_hh
