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

  /// Rebuild preconditioner every N Newton steps (1 = every step, default)
  CFuint m_lagFrequency;

  /// Rebuild when KSP iterations exceed this factor times the baseline count
  /// (0 = disabled, default 2.0)
  CFreal m_lagKSPGrowthThreshold;

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
