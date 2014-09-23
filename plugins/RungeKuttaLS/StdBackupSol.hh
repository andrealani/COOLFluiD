// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKuttaLS_StdBackupSol_hh
#define COOLFluiD_Numerics_RungeKuttaLS_StdBackupSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKLSData.hh"
#include "MathTools/RealVector.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This commands backs up the solution at the previous timestep
   * @author Kris Van den Abeele
   */
class StdBackupSol : public RKLSCom {
public:

  /**
   * Constructor.
   */
  explicit StdBackupSol(const std::string& name);

  /**
   * Destructor.
   */
  ~StdBackupSol()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for backup solution
  Framework::DataSocketSink<RealVector> socket_u0;

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class StdBackupSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKuttaLS_StdBackupSol_hh
