// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKutta2_StdSetup_hh
#define COOLFluiD_Numerics_RungeKutta2_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RK2Data.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "MathTools/RealVector.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to setup the RungeKutta2 Method
class StdSetup : public RK2Com {
public:

  /// Constructor.
  explicit StdSetup(const std::string& name);

  /// Destructor.
  ~StdSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for Rhs
  Framework::DataSocketSource<CFreal> socket_rhs;

  /// socket for backup solution
  Framework::DataSocketSource<RealVector> socket_u0;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSource<CFreal> socket_updateCoeff;

  /// socket for updateCoeff0
  /// back up of the first computed update coefficient
  Framework::DataSocketSource<CFreal> socket_upCoeff0;

  // handle to states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKutta2_StdSetup_hh

