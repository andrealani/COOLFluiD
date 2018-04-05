// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_BackwardEuler_StdSetup_hh
#define COOLFluiD_Numerics_BackwardEuler_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "BwdEulerData.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to setup the MeshData
/// according to the Method BwdEuler
/// @author Andrea Lani
/// @author Tiago Quintino
class StdSetup : public BwdEulerCom {
public:

  /// Constructor.
  explicit StdSetup(std::string name);

  /// Destructor.
  virtual ~StdSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  virtual void execute();

protected:

  /// socket for Rhs
  Framework::DataSocketSource<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSource<CFreal> socket_updateCoeff;

  /// socket for bStatesNeighbors
  /// It is a list of the neighbor states for the boundary states.
  /// It will be useful to avoid very expensive jacobian matrix
  /// reallocations when applying strong boundary condition
  Framework::DataSocketSource<std::valarray<Framework::State*> > socket_bStatesNeighbors;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BackwardEuler_StdSetup_hh

