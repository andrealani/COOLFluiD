// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKutta_StdBackupSol_hh
#define COOLFluiD_Numerics_RungeKutta_StdBackupSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKData.hh"
#include "MathTools/RealVector.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class StdBackupSol : public RKCom {
public:

  /// Constructor.
  explicit StdBackupSol(const std::string& name) :
    RKCom(name),
    socket_rhs("rhs"),
    socket_u0("u0"),
    socket_tempu("tempu"),
    socket_states("states")
  {
  }

  /// Destructor.
  ~StdBackupSol()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

  /// socket for backup solution
  Framework::DataSocketSink<
                              RealVector> socket_u0;

  /// socket for temporary storage of u
  Framework::DataSocketSink<
                              RealVector> socket_tempu;

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;


}; // class StdBackupSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKutta_StdBackupSol_hh

