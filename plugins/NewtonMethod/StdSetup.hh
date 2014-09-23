// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_StdSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// sent to Domain to be executed in order to setup the MeshData.
class StdSetup : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit StdSetup(const std::string& name);

  /// Destructor.
  ~StdSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for Rhs
  Framework::DataSocketSource<
                              CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSource<
                              CFreal> socket_updateCoeff;

  /// socket for PastStates's
  Framework::DataSocketSource<
                              Framework::State*> socket_pastStates;

  /// socket for bStatesNeighbors
  /// It is a list of the neighbor states for the boundary states.
  /// It will be useful to avoid very expensive jacobian matrix
  /// reallocations when applying strong boundary condition
  Framework::DataSocketSource<
                              std::valarray<Framework::State*> > socket_bStatesNeighbors;

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_StdSetup_hh
