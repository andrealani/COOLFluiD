// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ParJFSolveSys_hh
#define COOLFluiD_Numerics_Petsc_ParJFSolveSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/StdParSolveSys.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a standard command to be executed during
  * solution of the linear system using Petsc library
  *
  * @author Andrea Lani
  * @author Jiri Simonek
  *
  */
class ParJFSolveSys : public StdParSolveSys {
public:

  /**
   * Defines configuration options for preconditioner lagging.
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit ParJFSolveSys(const std::string& name);

  /**
   * Destructor.
   */
  ~ParJFSolveSys();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// socket for update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// How often to rebuild the preconditioner, counted in Newton iterations
  /// (NOT time steps): 1 (default) rebuilds at every iteration, i.e. no
  /// lagging; N rebuilds every N-th iteration; 0 disables this schedule, so
  /// rebuilds then happen only at the start of a new time step (see
  /// m_lagAcrossTimeSteps) or when the KSP count grows too much.
  /// Full explanation in defineConfigOptions() in the .cxx.
  CFuint m_lagFrequency;

  /// Safety net for lagging: rebuild as soon as a solve needs more than
  /// (this factor) x (the KSP iteration count right after the last rebuild),
  /// meaning the lagged preconditioner has degraded. 0 disables the check.
  /// Default 2.0.
  CFreal m_lagKSPGrowthThreshold;

  /// false (default): force a rebuild whenever a new (pseudo-)time step
  /// starts, the right choice for unsteady/time-accurate runs.
  /// true: let a lagged preconditioner survive across time steps. Needed
  /// for lagging in steady/pseudo-steady runs, where every Newton iteration
  /// counts as a new "time step" and the forced rebuild would otherwise
  /// override the LagFrequency schedule.
  bool m_lagAcrossTimeSteps;

  /// Call counter: increments every execute() call
  CFuint m_lagCounter;

  /// KSP iteration count from the most recent solve
  CFuint m_lastKSPIters;

  /// KSP iteration count at the time of the last preconditioner rebuild
  CFuint m_lastRebuildKSP;

  /// SubSystem iteration index at the time of the last preconditioner rebuild
  CFuint m_lastRebuildTimeStep;

}; // class SolveSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ParJFSolveSys_hh
