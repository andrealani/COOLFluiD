// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKuttaLS_RungeKuttaStep_hh
#define COOLFluiD_Numerics_RungeKuttaLS_RungeKuttaStep_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKLSData.hh"
#include "Framework/DataSocketSink.hh"
#include "MathTools/RealVector.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This command implements a R-K stage
   * @author Kris Van den Abeele
   */
class RungeKuttaStep : public RKLSCom {
public:

  /**
   * Constructor.
   */
  explicit RungeKuttaStep(const std::string& name);

  /**
   * Destructor.
   */
  ~RungeKuttaStep();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Setup private data of this class
   */
  virtual void setup();

  /**
   * Unsetup the private data of this class
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for backup solution
  Framework::DataSocketSink<RealVector> socket_u0;

  /// socket for updateCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKuttaLS_RungeKuttaStep_hh

