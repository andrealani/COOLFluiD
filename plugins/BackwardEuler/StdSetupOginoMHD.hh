// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_BackwardEuler_StdSetupOginoMHD_hh
#define COOLFluiD_Numerics_BackwardEuler_StdSetupOginoMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdSetup.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "MHD/MHDTerm.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Physics::MHD;

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class StdSetupOginoMHD : public StdSetup {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit StdSetupOginoMHD(const std::string& name);

  /// Destructor.
  ~StdSetupOginoMHD()
  {
  }

  /// Set up the private data
  void setup();

  /// Execute Processing actions
  void execute();

  /// Configures this object by complementing the
  /// implementation in ConfigObject
  void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private:

  /// socket for flags telling if the corresponding state has to be updated
  Framework::DataSocketSource< bool> socket_hasToBeUpdated;

  /// vector of state coordinates
  RealVector _stateCoord;

  /// radius of the planet
  CFreal _radiusPlanet;

  /// radius of the transition region which is the limit of effect of the planetary state
  /// values
  CFreal _radiusTransition;

  /// the solution variables vector inside the planet
  std::vector< CFreal > _planetSolutionVars;

}; // class StdSetupOginoMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BackwardEuler_StdSetupOginoMHD_hh
