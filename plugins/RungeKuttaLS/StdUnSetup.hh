// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKuttaLS_StdUnSetup_hh
#define COOLFluiD_Numerics_RungeKuttaLS_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKLSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a unsets the RungeKuttaLS method.
 *
 * @author Kris Van den Abeele
 */
class StdUnSetup : public RKLSCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdUnSetup()
  {
  }

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected:

  /// socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for backup solution
  Framework::DataSocketSink<RealVector> socket_u0;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKuttaLS_StdUnSetup_hh

